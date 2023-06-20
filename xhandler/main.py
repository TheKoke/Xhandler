import os, sys
import threading

from analysis import Analytics, PeakSupervisor
from workbook import WorkbookMaster, WorkbookParser
from shunting_yard import ReactionMaster, Reaction

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
            SELECTED_DOTS_X.append(event.xdata)
            SELECTED_DOTS_Y.append(event.ydata)

            print(f'Point (x: {round(event.xdata, 2)}, y: {round(event.ydata, 2)}) was selected.')
        
    def clear_all_stuff(self) -> None:
        self.clear_selected_dots()
        self.axes.clear()
        self.figure.canvas.draw()

    def clear_selected_dots(self) -> None:
        SELECTED_DOTS_X.clear()
        SELECTED_DOTS_Y.clear()

    def draw_uncalibrated_spectrum(self, spectrum: np.ndarray) -> None:
        self.axes.clear()

        self.axes.plot(np.arange(1, len(spectrum) + 1), spectrum)
        self.figure.canvas.draw() 

    def draw_calibrated_spectrum(self, spectrum: np.ndarray, energy_view: np.ndarray) -> None:
        self.axes.clear()

        self.axes.plot(energy_view, spectrum)
        self.figure.canvas.draw()

    def draw_peak(self, peak: PeakSupervisor) -> None:
        curve = peak.lorentzian

        self.axes.plot(curve.three_sigma(), curve.pdf(), color='red')
        self.figure.canvas.draw()

    def scat_dots(self, xs: np.ndarray, ys: np.ndarray) -> None:
        self.axes.scatter(xs, ys, color='red', label='Theoretical peaks center.')
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
            msg = self.handle_command(comm)
            if msg == 'quit':
                break

            print(msg)

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
            case 'quit': return 'quit'
            case _: return self.error_message()

    def take_reaction(self) -> Reaction:
        print('We need info about your nuclear reaction.')
        str_react = input('Please write down analyzing nuclear reaction: ')
        energy = float(input('Type here beam energy (in MeV): '))

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
        prepared = next(i for i in CACHED_SPECTRES if i.angle == angle)
        if not prepared.is_calibrated:
            self.observer.draw_uncalibrated_spectrum(prepared.spectrum)

        self.observer.draw_calibrated_spectrum(prepared.spectrum, prepared.energy_view())
        for p in prepared.peaks:
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

        if answer.isdigit() and float(answer) in [i.angle for i in CACHED_SPECTRES]:
            analyzed = next(i for i in CACHED_SPECTRES if i.angle == float(answer))
            self.workbooker.write(str(analyzed) + '\n\n')
            return 'Analyzed parameters was wroted to workbook.'
        else:
            return 'Cannot find this angle inside the analyzed ones.'

    def __show_analyzed_spectres(self) -> None:
        cached = [i.angle for i in CACHED_SPECTRES]

        print('Choose the angle in analyzed spectrums angles:')
        print(cached)

    def change(self) -> str:
        pass

    def read(self) -> str:
        pass

    def calibrate(self) -> str:
        if not self.is_spectrum_opened:
            return 'First open the spectrum.'

        print('You can double click to graph, to pick points on them.')
        print('The last picked 2 points will be used for calibration.')
        input('Please, press enter to continue, when you finish selecting of points.\n')

        val, e0 = self.analitics.calibrate((int(SELECTED_DOTS_X[-1]), int(SELECTED_DOTS_X[-2])))

        self.observer.draw_calibrated_spectrum(self.analitics.spectrum, self.analitics.energy_view())

        peaks_indexes = self.analitics.try_find_peaks()
        self.observer.scat_dots(self.analitics.theory_peaks[:len(peaks_indexes)], self.analitics.spectrum[peaks_indexes])

        return 'calibrated by: ' + \
        f'E(ch) = {round(val, 3)}ch + {round(e0, 3)}'

    def fit_peak(self) -> str:
        if not self.is_spectrum_opened:
            return 'First open the spectrum.'

        if not self.analitics.is_calibrated:
            return 'Spectrum must be calibrated first.'
        
        peaks = self.analitics.create_peaks()
        for peak in peaks:
            self.observer.draw_peak(peak)

        return 'All peaks was drown.'

    def save(self) -> str:
        CACHED_SPECTRES.append(self.analitics)

        self.analitics = None
        self.is_spectrum_opened = False

        return f'Spectrum of {CACHED_SPECTRES[-1].angle} degree was saved.\n' + \
                'To write this to workbook type *write down*'

    def error_message(self) -> str:
        return '404 Error... Command not found.'


def startup() -> None:
    print('!!!!!START PROCESSING!!!!!')
    print("Welcome to command-promt software to analyze nuclear reaction's spectrums.\n\n")

    path = input('Please write the PATH to spectres: ')
    main_commandor = Commandor(path)

    loop = threading.Thread(target=main_commandor.main)
    loop.start()
    pyplot.show()


if __name__ == '__main__':
    startup()
