import numpy as np
from maths import *
from physics import *


class Peak:
    def __init__(self, x_data: np.ndarray, y_data: np.ndarray, center_index: int) -> None:
        self.x_data = x_data
        self.y_data = y_data

        self.mu_index = center_index
        self.mu = self.x_data[center_index]

        self.PEAK_WIDTH = len(x_data) // 50

        self.gaussian = self.approximate()

    def approximate(self) -> Gaussian:
        return Gaussian(
            np.vstack(
                (self.x_data[self.mu_index - self.PEAK_WIDTH // 2: self.mu_index + self.PEAK_WIDTH // 2],
                self.y_data[self.mu_index - self.PEAK_WIDTH // 2: self.mu_index + self.PEAK_WIDTH // 2])
            )
        )


class Analytics:
    def __init__(self, spectrum: np.ndarray, reaction: Reaction, angle: float) -> None:
        self.angle = angle
        self.spectrum = spectrum

        self.reaction = reaction
        self.states = reaction.residual.states

        self.theory_peaks = [self.reaction.fragment_energy(state, angle) for state in self.states]

        self.is_calibrated = False
        self.scale_value, self.scale_shift = 0, 0
        self.peaks: list[Peak] = []

    def __str__(self) -> str:
        self.peaks = sorted(self.peaks, key=lambda peak: peak.mu, reverse=True)

        info = f'{self.angle} - angle spectrum analysis: \n'
        info += 'Calibrated by equation: E(ch) = ' + \
        f'{round(self.scale_value, 3)} * ch + {round(self.scale_shift, 3)}\n'

        info += f'--Peaks analysis info--'.center(66) + '\n'
        info += 'Carbon state, MeV'.center(20) + '\t' + 'center, MeV'.center(15) + '\t' 
        info += 'fwhm, MeV'.center(15) + '\t' + 'area'.center(15) + '\n'
        for i in range(len(self.peaks)):
            info += str(round(self.states[i], 3)).center(20) + '\t'
            info += str(round(self.peaks[i].mu, 3)).center(15) + '\t'
            info += str(round(self.peaks[i].gaussian.fwhm(), 3)).center(15) + '\t'
            info += str(round(self.peaks[i].gaussian.area, 3)).center(15) + '\n'

        return info

    def calibrate(self, anchors: tuple[int, int]) -> tuple[float, float]:
        anchors = sorted(anchors, reverse=True)

        matrix = np.array([[anchors[0], 1], [anchors[1], 1]])
        right_side = np.array([self.theory_peaks[0], self.theory_peaks[1]])
        solution = np.linalg.solve(matrix, right_side)

        self.scale_value, self.scale_shift = solution[0], solution[1]
        self.is_calibrated = True

        return (self.scale_value, self.scale_shift)
    
    def try_find_peaks(self) -> list[int]:
        if not self.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before finding peaks.')

        return [int((theory - self.scale_shift) / self.scale_value) for theory in self.theory_peaks]

    def create_peak(self, center: float) -> Peak:
        if not self.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before creating peaks.')

        center_index = int((center - self.scale_shift) / self.scale_value)
        energy_view = np.arange(1, len(self.spectrum) + 1) * self.scale_value + self.scale_shift

        self.peaks.append(Peak(energy_view, self.spectrum, center_index))
        return self.peaks[-1]


if __name__ == '__main__':
    pass
