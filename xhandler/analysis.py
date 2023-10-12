import numpy as np
from base import Spectrum, Peak
from maths import Gaussian
from physics import Reaction


class PeakSupervisor:
    def __init__(self, x_data: np.ndarray, y_data: np.ndarray, center: float, fwhm: float) -> None:
        self.x_data = x_data
        self.y_data = y_data

        self.mu = center
        self.mu_index = np.abs(np.subtract(x_data, center)).argmin()

        self.fwhm = fwhm
        self.peak = Peak()
        self.lorentzian = self.approximate()

    def approximate(self) -> Gaussian:
        peak_start = self.mu_index - self.width() // 2
        peak_stop = self.mu_index + self.width() // 2

        peak_start = peak_start if peak_start >= 0 else 0
        peak_stop - peak_stop if peak_stop < len(self.y_data) else len(self.y_data) - 1

        area = self.tuck_up_area(self.x_data[peak_start: peak_stop], self.y_data[peak_start: peak_stop])
        self.__save_params(self.mu, self.fwhm, area)
        return Gaussian(self.mu, self.fwhm, area)
    
    def width(self) -> int:
        scale_value = (self.x_data[-1] - self.x_data[0]) / len(self.x_data)
        tenth_width = self.fwhm / np.log10(2)
        return int(tenth_width / scale_value)
    
    def tuck_up_area(self, xdata: np.ndarray, ydata: np.ndarray) -> float:
        # xs = 2 / np.pi * (self.fwhm / (4 * (xdata - self.mu) ** 2 + self.fwhm ** 2))
        sigma = self.fwhm / (2 * np.sqrt(2 * np.log(2)))
        xs = 1 / np.sqrt(2 * np.pi * sigma ** 2) * np.exp(- (xdata - self.mu) ** 2 / (2 * sigma ** 2))
        return (ydata * xs).sum() / (xs ** 2).sum()
    
    def __save_params(self, mu: float, fwhm: float, area: float) -> None:
        self.peak.mu = mu
        self.peak.fwhm = fwhm
        self.peak.area = area


class Analytics:
    def __init__(self, spectrum: np.ndarray, reaction: Reaction, angle: float) -> None:
        self.base = Spectrum()

        self.base.angle = angle
        self.base.spectrum = spectrum

        self.reaction = reaction
        self.theory_peaks = self.found_theory_peaks()

        self.peaks: list[PeakSupervisor] = []
        self.truncate_spectrum()

    @property
    def angle(self) -> float:
        return self.base.angle

    @property
    def spectrum(self) -> np.ndarray:
        return self.base.spectrum

    @property
    def is_calibrated(self) -> bool:
        return self.base.is_calibrated

    def __str__(self) -> str:
        self.peaks = sorted(self.peaks, key=lambda peak: peak.mu, reverse=True)

        info = f'{self.base.angle} - angle spectrum analysis: \n'
        info += 'Calibrated by equation: E(ch) = ' + \
        f'{round(self.base.scale_value, 3)} * ch + {round(self.base.scale_shift, 3)}\n'

        info += f'--Peaks analysis info--'.center(66) + '\n'
        info += 'Fragment state, MeV'.center(20) + '\t' + 'center, MeV'.center(15) + '\t' 
        info += 'fwhm, MeV'.center(15) + '\t' + 'area'.center(15) + '\n'
        for i in range(len(self.peaks)):
            info += str(round(self.reaction.residual.states[i], 3)).center(20) + '\t'
            info += str(round(self.peaks[i].mu, 3)).center(15) + '\t'
            info += str(round(self.peaks[i].fwhm, 3)).center(15) + '\t'
            info += str(round(self.peaks[i].peak.area, 3)).center(15) + '\n'

        return info
    
    def calibrate(self, anchors: tuple[int, int]) -> tuple[float, float]:
        anchors = sorted(anchors, reverse=True)

        matrix = np.array([[anchors[0], 1], [anchors[1], 1]])
        right_side = np.array([self.theory_peaks[0], self.theory_peaks[1]])
        solution = np.linalg.solve(matrix, right_side)

        self.base.scale_value, self.base.scale_shift = solution[0], solution[1]
        return (self.base.scale_value, self.base.scale_shift)
    
    def energy_view(self) -> np.ndarray:
        return self.base.energy_view
    
    # TODO: Fix this crouch. Needs little bit architectural touch.
    def found_theory_peaks(self) -> list[float]:
        collected = []
        for state in self.reaction.residual.states:
            energy_after_reaction = self.reaction.fragment_energy(state, self.base.angle)
            collected.append(energy_after_reaction)

        return collected
    
    def try_find_peaks(self) -> list[int]:
        if not self.base.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before finding peaks.')
        
        visible = []
        for theory in self.theory_peaks:
            pretend = int((theory - self.base.scale_shift) / self.base.scale_value) - 1
            if pretend <= 0:
                continue

            visible.append(pretend)

        return visible

    def create_peaks(self) -> list[PeakSupervisor]:
        if not self.base.is_calibrated:
            raise RuntimeError('Spectrum must be calibrated before creating peaks.')

        theory_indexes = self.try_find_peaks()
        for i in range(len(theory_indexes)):
            xs = self.base.energy_view
            ys = self.base.spectrum
            center = self.energy_view()[theory_indexes[i]]
            gamma = self.reaction.residual.wigner_widths[i]

            current = PeakSupervisor(xs, ys, center, gamma)
            self.peaks.append(current)

        self.base.peaks = [supervisor.peak for supervisor in self.peaks]
        return self.peaks
    
    def truncate_spectrum(self) -> None:
        count = 0
        for i in range(len(self.base.spectrum)):
            if count >= 50:
                self.base.spectrum = self.base.spectrum[:i - 20]
                break

            if self.base.spectrum[i] == 0:
                count += 1
            else:
                count = 0


class Sectioner:
    def __init__(self, reaction: Reaction, spectres: list[Spectrum]) -> None:
        self.reaction = reaction
        self.spectres = spectres

    def centermass_cs(self) -> np.ndarray:
        pass

    def lab_cs(self) -> np.ndarray:
        pass

    def centermass_angles(self) -> np.ndarray:
        pass

    def __x_square(self) -> np.ndarray:
        pass

    def __g_constant(self) -> float:
        pass


if __name__ == '__main__':
    pass
