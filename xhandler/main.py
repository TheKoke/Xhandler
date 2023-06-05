import os
import threading

from analysis import *
from workbook import *
from shunting_yard import *

import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib.backend_bases import MouseEvent


SELECTED_DOTS_X = []
SELECTED_DOTS_Y = []
CACHED_SPECTRES: list[Analytics] = []


class Sleuth:
    def __init__(self, directory: str) -> None:
        self.dir = directory
        self.angles = self.get_angles_list()

    def get_angles_list(self) -> list[float]:
        files = self.only_files(os.listdir(self.dir))
        return [self.angle_from_name(file) for file in files]

    def angle_from_name(self, file: str) -> float:
        return float(file.split('.')[0])
    
    def only_files(self, dirs: list[str]) -> list[str]:
        return [direc for direc in dirs if '.' in direc]
    
    def to_spectrum(self, angle: float) -> np.ndarray:
        return np.loadtxt(f'{self.dir}\\{int(angle)}.txt')


class Observer:
    def __init__(self) -> None:
        self.figure, self.axes = pyplot.subplots()
        self.figure.canvas.mpl_connect('button_press_event', self.select_point)

    def select_point(self, event: MouseEvent) -> str | None:
        if event.dblclick:
            if len(SELECTED_DOTS_X) > 10: SELECTED_DOTS_X.clear()
            if len(SELECTED_DOTS_Y) > 10: SELECTED_DOTS_Y.clear()
            
            SELECTED_DOTS_X.append(event.xdata)
            SELECTED_DOTS_Y.append(event.ydata)

            print(f'Point (x: {round(event.xdata, 2)}, y: {round(event.ydata, 2)}) was selected.')
        
    def clear_all_stuff(self) -> None:
        self.axes.clear()
        self.figure.canvas.draw()

    def draw_uncalibrated_spectrum(self, spectrum: np.ndarray) -> None:
        self.axes.clear()

        self.axes.plot(np.arange(1, len(spectrum) + 1), spectrum)
        self.figure.canvas.draw() 

    def draw_calibrated_spectrum(self, spectrum: np.ndarray, energy_view: np.ndarray) -> None:
        self.axes.clear()

        self.axes.plot(energy_view, spectrum)
        self.figure.canvas.draw()

    def draw_peak(self, peak: Peak) -> None:
        gauss = peak.gaussian

        self.axes.plot(gauss.three_sigma(), gauss.function(), color='red')
        self.figure.canvas.draw()

#TODO: Implement all class.
class Commandor:
    def __init__(self, path_to_find: str) -> None:
        self.sleuth = Sleuth(path_to_find)
        self.observer = Observer()
        self.workbooker = WorkbookMaster(path_to_find)

        self.analitics: Analytics = None
        self.reaction: Reaction = None

        self.is_spectrum_opened = False

    def main(self) -> None:
        self.reaction = self.take_reaction()
        while True:
            print(self.show_possible_commands())

            comm = input('What do we do? Type here: ')
            print(self.handle_command(comm))

    def show_possible_commands(self) -> list[str]:
        if not self.is_spectrum_opened:
            return ['open', 'change', 'read', 'write down', 'quit']
        else:
            return ['calibrate', 'fit peak', 'save', 'close']

    def handle_command(self, command: str) -> str:
        match command:
            case 'open': return self.open()
            case 'close': return self.close()
            case 'write down': return self.write_down()
            case 'change': return self.change()
            case 'read': return self.read()
            case 'calibrate': return self.calibrate()
            case 'fit peak': return self.fit_peak()
            case 'save': return self.save()
            case _: return self.error_message()

    def take_reaction(self) -> Reaction:
        print('We need info about your nuclear reaction.')
        str_react = input('Please write down analyzing nuclear reaction: ')
        energy = float(input('Type here beam energy: '))

        return ReactionMaster(str_react, energy).to_reaction()

    def open(self) -> str:
        print('Finded angles:')
        print(self.sleuth.angles)

        comm = input('Choose the angle to show they spectrum: ')
        if comm == 'cancel':
            return 'Opening operation was cancelled.'
        
        if not comm.isdigit():
            print('Please, write correctly!')
            return self.open()
        
        self.__open_spectrum(comm)
        return f'Spectrum of {self.analitics.angle} was opened.'
    
    def __open_spectrum(self, pretend: str) -> None:
        pretend = float(pretend)
        if pretend in CACHED_SPECTRES:
            self.__open_cached_spectrum(pretend)
            return

        spectrum = self.sleuth.to_spectrum(pretend)

        self.observer.draw_uncalibrated_spectrum(spectrum)

        self.analitics = Analytics(spectrum, self.reaction, pretend)
        self.is_spectrum_opened = True

    def __open_cached_spectrum(self, angle: float) -> None:
        analysis = next(i for i in CACHED_SPECTRES if i.angle == angle)
        if not analysis.is_calibrated:
            self.observer.draw_uncalibrated_spectrum(analysis.spectrum)

        energy_view = analysis.scale_value * np.arange(1, len(analysis.spectrum) + 1)
        energy_view += analysis.scale_shift

        self.observer.draw_calibrated_spectrum(analysis.spectrum, energy_view)
        for p in analysis.peaks:
            self.observer.draw_peak(p)

    def close(self) -> str:
        print("If you quit immediately, your changes doesn't applies.")
        answer = input('Are you serious? Yes or No: ')

        if answer.lower() == 'yes':
            self.observer.clear_all_stuff()
            self.analitics = None
            self.is_spectrum_opened = False

            return f'The spectrum was closed with no save.'
        else:
            return 'That is a big deal. Save it to apply your changes.'

    def write_down(self) -> str:
        self.__show_analyzed_spectres()

        answer = input('Type here: ')

        if answer.isdigit() and answer in self.workbooker.gather_analyzed():
            self.workbooker.write(str(self.analitics))
            return 'Analyzed parameters was wroted to workbook.'
        else:
            return 'Cannot find this angle inside the analyzed ones.'

    def __show_analyzed_spectres(self) -> None:
        print('Choose the angle in analyzed spectrums angles:')
        print(self.workbooker.gather_analyzed())

    def change(self) -> str:
        pass

    def read(self) -> str:
        pass

    def calibrate(self) -> str:
        print('You can double click to graph, to pick points on them.')
        print('The last picked 2 points will be used for calibration.')
        input('Please, press enter to continue, when you finish selecting of points.\n')

        self.analitics.calibrate((int(SELECTED_DOTS_X[-1]), int(SELECTED_DOTS_X[-2])))
        self.__clear_selected_dots()

        energy_view = self.analitics.scale_value * np.arange(1, len(self.analitics.spectrum) + 1)
        energy_view += self.analitics.scale_shift
        self.observer.draw_calibrated_spectrum(self.analitics.spectrum, energy_view)

        return 'calibrated by: ' + \
        f'E(ch) = {round(self.analitics.scale_value, 3)}ch + {round(self.analitics.scale_shift, 3)}'

    def fit_peak(self) -> str:
        if not self.analitics.is_calibrated:
            return 'Spectrum must be calibrated first.'

        print('Double click for pick center of peak.')
        print('This software can analyze only one peak in a row.')
        input('Please, press the enter to continue.\n')

        created = self.analitics.create_peak(SELECTED_DOTS_X[-1])
        self.observer.draw_peak(created)

        return f'Peak at {created.mu} MeV was drawed.'

    def save(self) -> str:
        CACHED_SPECTRES.append(self.analitics)

        self.analitics = None
        self.is_spectrum_opened = False

        return f'Spectrum of {CACHED_SPECTRES[-1].angle} degrees was saved.' + \
                'To write this to workbook type *write down*'

    def error_message(self) -> str:
        return '404 Error... Command not found.'
    
    def __clear_selected_dots(self) -> None:
        SELECTED_DOTS_X.clear()
        SELECTED_DOTS_Y.clear()


def startup() -> None:
    print('!!!!!START PROCESSING!!!!!')
    print("Welcome to command-promt software to analyze nuclear reaction's spectrums.\n\n")
    # print('This version for C12 + B10 reaction, for more download the full-version.\n\n')

    path = input('Please write the PATH to spectres: ')
    main_commandor = Commandor(path)

    loop = threading.Thread(target=main_commandor.main)
    loop.start()
    pyplot.show()


if __name__ == '__main__':
    startup()
