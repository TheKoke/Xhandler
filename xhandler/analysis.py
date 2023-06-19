import numpy as np
from base import Notebook, Peak
from maths import Lorentzian
from physics import Reaction, Nuclei, Ionization


class PeakSupervisor:
    def __init__(self, x_data: np.ndarray, y_data: np.ndarray, center: float, fwhm: float) -> None:
        self.x_data = x_data
        self.y_data = y_data

        self.mu = center
        self.mu_index = np.abs(np.subtract(x_data, center)).argmin()

        self.fwhm = fwhm
        self.lorentzian = self.approximate()

        self.base = Peak()

    def approximate(self) -> Lorentzian:
        peak_start = self.mu_index - self.width() // 2
        peak_stop = self.mu_index + self.width() // 2

        peak_start = peak_start if peak_start >= 0 else 0
        peak_stop - peak_stop if peak_stop < len(self.y_data) else len(self.y_data) - 1

        area = self.tuck_up_area(self.x_data[peak_start: peak_stop], self.y_data[peak_start: peak_stop])
        self.__save_params(self.mu, self.fwhm, area)
        return Lorentzian(self.mu, self.fwhm, area)
    
    def width(self) -> int:
        scale_value = (self.x_data[-1] - self.x_data[0]) / len(self.x_data)
        tenth_width = self.fwhm / np.log10(2)
        return int(tenth_width / scale_value)
    
    def tuck_up_area(self, xdata: np.ndarray, ydata: np.ndarray) -> float:
        xs = 2 / np.pi * (self.fwhm / (4 * (xdata - self.mu) ** 2 + self.fwhm ** 2))
        return (ydata * xs).sum() / (xs ** 2).sum()
    
    def __save_params(self, mu: float, fwhm: float, area: float) -> None:
        self.base.mu = mu
        self.base.fwhm = fwhm
        self.base.area = area


class Analytics:
    def __init__(self, spectrum: np.ndarray, reaction: Reaction, angle: float) -> None:
        self.basenote = Notebook()

        self.basenote.angle = angle
        self.basenote.spectrum = spectrum

        self.reaction = reaction
        self.theory_peaks = self.found_theory_peaks()

        self.peaks: list[Peak] = []
        self.truncate_spectrum()

    def __str__(self) -> str:
        self.peaks = sorted(self.peaks, key=lambda peak: peak.mu, reverse=True)

        info = f'{self.basenote.spectrum} - angle spectrum analysis: \n'
        info += 'Calibrated by equation: E(ch) = ' + \
        f'{round(self.basenote.scale_value, 3)} * ch + {round(self.basenote.scale_shift, 3)}\n'

        info += f'--Peaks analysis info--'.center(66) + '\n'
        info += 'Fragment state, MeV'.center(20) + '\t' + 'center, MeV'.center(15) + '\t' 
        info += 'fwhm, MeV'.center(15) + '\t' + 'area'.center(15) + '\n'
        for i in range(len(self.peaks)):
            info += str(round(self.reaction.residual.states[i], 3)).center(20) + '\t'
            info += str(round(self.peaks[i].mu, 3)).center(15) + '\t'
            info += str(round(self.peaks[i].fwhm, 3)).center(15) + '\t'
            info += str(round(self.peaks[i].area, 3)).center(15) + '\n'

        return info
    
    def calibrate(self, anchors: tuple[int, int]) -> tuple[float, float]:
        anchors = sorted(anchors, reverse=True)

        matrix = np.array([[anchors[0], 1], [anchors[1], 1]])
        right_side = np.array([self.theory_peaks[0], self.theory_peaks[1]])
        solution = np.linalg.solve(matrix, right_side)

        self.basenote.scale_value, self.basenote.scale_shift = solution[0], solution[1]
        return (self.basenote.scale_value, self.basenote.scale_shift)
    
    # TODO: Fix this crouch. Needs little bit architectural touch.
    def found_theory_peaks(self) -> list[float]:
        bete_bloch = Ionization(self.reaction.fragment, Nuclei(58, 34), detector='C4H10')

        collected = []
        for state in self.reaction.residual.states:
            energy_after_reaction = self.reaction.fragment_energy(state, self.basenote.angle)
            energy_loss = bete_bloch.energy_loss(energy_after_reaction, 4)

            collected.append(energy_after_reaction - energy_loss)

        return collected
    
    def try_find_peaks(self) -> list[int]:
        if not self.basenote.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before finding peaks.')
        
        visible = []
        for theory in self.theory_peaks:
            pretend = int((theory - self.basenote.scale_shift) / self.basenote.scale_value) - 1
            if pretend < 0:
                continue

            visible.append(pretend)

        return visible

    def create_peaks(self) -> list[Peak]:
        if not self.basenote.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before creating peaks.')

        for i in range(len(self.theory_peaks)):
            xs = self.basenote.energy_view
            ys = self.basenote.spectrum
            center = self.theory_peaks[i]
            gamma = self.reaction.residual.wigner_width[i]

            current = PeakSupervisor(xs, ys, center, gamma)
            self.peaks.append(current.base)

        self.basenote.peaks = self.peaks.copy()
        return self.peaks.copy()
    
    def truncate_spectrum(self) -> None:
        count = 0
        for i in range(len(self.basenote.spectrum)):
            if count >= 50:
                self.basenote.spectrum = self.basenote.spectrum[:i - 20]
                break

            if self.basenote.spectrum[i] == 0:
                count += 1
            else:
                count = 0


if __name__ == '__main__':
    pass
