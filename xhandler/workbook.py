import os
import numpy as np

from analysis import *


class WorkbookParser:
    def __init__(self, path: str, spectres_path: str) -> None:
        self.path = path
        self.spectres_path = spectres_path

    def collect_all_prepared(self) -> list[Analytics]:
        all_reports = open(self.spectres_path, 'r').read().split('\n\n')
        return [self.reverse_parameters(report) for report in all_reports]

    def reverse_parameters(self, report: str) -> Analytics:
        result = Analytics(self.__get_spectrum(self.spectres_path, self.__get_angle(report)))

        result.scale_value, result.scale_shift = self.__get_calibration_constants(report)
        result.peaks = self.__get_peaks(report)
        result.is_calibrated = True

        return result

    def __get_spectrum(self, path: str, angle: float) -> np.ndarray:
        return np.loadtxt(path + f'{int(angle)}.txt')

    def __get_angle(self, report: str) -> float:
        pass

    def __get_calibration_constants(self, report: str) -> tuple[float, float]:
        pass

    def __get_peaks(self, report: str) -> list[Peak]:
        pass


class WorkbookMaster:
    def __init__(self, spectres_path: str) -> None:
        self.path = os.getcwd() + '\\workbook.txt'
        self.spectres_path = spectres_path
        self.analyzed_spectres: list[Analytics] = []


if __name__ == '__main__':
    pass
