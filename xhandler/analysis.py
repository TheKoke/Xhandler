import numpy as np
from maths import *
from physics import *


class Peak:
    def __init__(self, x_data: np.ndarray, y_data: np.ndarray, center_index: int, width: int) -> None:
        self.x_data = x_data
        self.y_data = y_data

        self.mu_index = center_index
        self.mu = self.x_data[center_index]

        self.width = width
        self.gaussian = self.approximate()

    def approximate(self) -> Gaussian:
        peak_start = self.mu_index - self.width // 2
        peak_stop = self.mu_index + self.width // 2

        peak_start = peak_start if peak_start >= 0 else 0
        peak_stop - peak_stop if peak_stop < len(self.y_data) else len(self.y_data) - 1

        return Gaussian(np.vstack((self.x_data[peak_start: peak_stop], self.y_data[peak_start: peak_stop])))


class Analytics:
    def __init__(self, spectrum: np.ndarray, reaction: Reaction, angle: float) -> None:
        self.angle = angle
        self.spectrum = spectrum

        self.reaction = reaction
        self.states = reaction.residual.states

        self.theory_peaks = self.found_theory_peaks()

        self.is_calibrated = False
        self.scale_value, self.scale_shift = 0, 0
        self.peaks: list[Peak] = []

        self.truncate_spectrum()

    def __str__(self) -> str:
        self.peaks = sorted(self.peaks, key=lambda peak: peak.mu, reverse=True)

        info = f'{self.angle} - angle spectrum analysis: \n'
        info += 'Calibrated by equation: E(ch) = ' + \
        f'{round(self.scale_value, 3)} * ch + {round(self.scale_shift, 3)}\n'

        info += f'--Peaks analysis info--'.center(66) + '\n'
        info += 'Fragment state, MeV'.center(20) + '\t' + 'center, MeV'.center(15) + '\t' 
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
    
    # TODO: Fix this crouch. Needs little bit architectural touch.
    def found_theory_peaks(self) -> list[float]:
        bete_bloch = Ionization(self.reaction.fragment, Nuclei(58, 34), detector='C4H10')

        collected = []
        for state in self.states:
            energy_after_reaction = self.reaction.fragment_energy(state, self.angle)
            energy_loss = bete_bloch.energy_loss(energy_after_reaction, 4)

            collected.append(energy_after_reaction - energy_loss)

        return collected
    
    def try_find_peaks(self) -> list[int]:
        if not self.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before finding peaks.')
        
        visible = []
        for theory in self.theory_peaks:
            pretend = int((theory - self.scale_shift) / self.scale_value) - 1
            if pretend < 0:
                continue

            visible.append(pretend)

        return visible

    def create_peaks(self) -> list[Peak]:
        if not self.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before creating peaks.')

        peak_width = self.define_peak_width()
        energy_view = self.scale_shift + self.scale_value * np.arange(1, len(self.spectrum) + 1)

        peaks = self.try_find_peaks()
        for peak in peaks:
            self.peaks.append(Peak(energy_view, self.spectrum, peak, peak_width))

        return self.peaks.copy()

    # TODO: Fix this croutch. Also needs architectural touch.
    def define_peak_width(self) -> int:
        fwhm = 0.826
        tenth_width = fwhm / np.log10(2)
        return int(tenth_width / self.scale_value)
    
    def truncate_spectrum(self) -> None:
        count = 0
        for i in range(len(self.spectrum)):
            if count >= 50:
                self.spectrum = self.spectrum[:i - 20]
                break

            if self.spectrum[i] == 0:
                count += 1
            else:
                count = 0


if __name__ == '__main__':
    pass
