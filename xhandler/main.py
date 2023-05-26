import os
import threading

from analysis import *
from workbook import *

import numpy as np
import matplotlib.pyplot as pyplot
from matplotlib.backend_bases import MouseEvent


SELECTED_DOTS_X = []
SELECTED_DOTS_Y = []


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

            return f'Point (x: {round(event.xdata, 2)}, y: {round(event.ydata, 2)}) was selected.'
        
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


class Commandor:
    def __init__(self, path_to_find: str) -> None:
        self.sleuth = Sleuth(path_to_find)
        self.observer = Observer()
        self.workbooker = WorkbookMaster(path_to_find + '\\workbook.txt')

        self.analitics: Analytics = None

        self.is_spectrum_opened = False

    def main(self) -> None:
        while True:
            print(self.show_possible_commands())

            comm = input('What do we do? Type here: ')
            print(self.handle_command(comm))

    def show_possible_commands(self) -> list[str]:
        if not self.is_spectrum_opened:
            return ['open', 'change', 'read', 'write down']
        else:
            return ['calibrate', 'fit peak', 'save']

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
        spectrum = self.sleuth.to_spectrum()

        self.observer.draw_uncalibrated_spectrum(spectrum)

        self.analitics = Analytics(spectrum, pretend)
        self.is_spectrum_opened = True

    def close(self) -> str:
        print("I fyou quit immediately, your changes doesn't applies.")
        answer = input('Are you serious? Yes or No: ')

        if answer.lower() == 'yes':
            self.observer.clear_all_stuff()
            return f'Spectrum of {self.analitics.angle} was closed with no save.'
        else:
            return 'That is a big deal. Save it to apply your changes.'

    def write_down(self) -> str:
        pass

    def change(self) -> str:
        pass

    def read(self) -> str:
        pass

    def calibrate(self) -> str:
        pass

    def fit_peak(self) -> str:
        pass

    def save(self) -> str:
        pass

    def error_message(self) -> str:
        return '404 Error... Command not found.'
    
    def __collect_analyzed(self) -> list[Analytics]:
        pass


def startup() -> None:
    print('!!!!!START PROCESSING!!!!!')
    print("Welcome to command-promt software to analyze nuclear reaction's spectrums.")
    print('This version for C12 + B10 reaction, for more download the full-version.\n\n')

    path = os.getcwd() + '\\fragment'
    main_commandor = Commandor(path)

    loop = threading.Thread(target=main_commandor.main)
    loop.start()
    pyplot.show()


if __name__ == '__main__':
    startup()
