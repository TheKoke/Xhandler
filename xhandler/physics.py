import numpy as np
from informer import Informator


class Nuclei:
    def __init__(self, charge: int, nuclons: int) -> None:
        self.nuclons = nuclons
        self.charge = charge

    def __str__(self) -> str:
        return Informator.name(self.charge, self.nuclons)

    @property
    def mass_excess(self) -> float:
        return Informator.mass_excess(self.charge, self.nuclons)
    
    @property
    def states(self) -> list[float]:
        return Informator.states(self.charge, self.nuclons)
    
    @property
    def wigner_widths(self) -> list[float]:
        return Informator.wigner_widths(self.charge, self.nuclons)
    
    @property
    def mass(self) -> float:
        return self.charge * 938.27 + (self.nuclons - self.charge) * 939.57


class Reaction:
    def __init__(self, beam: Nuclei, target: Nuclei, fragment: Nuclei, beam_energy: float) -> None:
        self.beam = beam
        self.target = target
        self.fragment = fragment
        self.residual = self.__residual_nuclei()

        self.beam_energy = beam_energy

    @property
    def is_elastic(self) -> bool:
        return self.beam == self.fragment

    def __residual_nuclei(self) -> Nuclei:
        nuclon = (self.beam.nuclons + self.target.nuclons) - self.fragment.nuclons
        charge = (self.beam.charge + self.target.charge) - self.fragment.charge
        return Nuclei(charge, nuclon)

    def reaction_quit(self, residual_state: float = 0) -> float:
        q0 = (self.beam.mass_excess + self.target.mass_excess) - (self.fragment.mass_excess + self.residual.mass_excess)
        return q0 - residual_state
    
    def fragment_energy(self, residual_state: float, fragment_angle: float) -> np.ndarray:
        r = Reaction.__r_factor(
            self.beam.mass, 
            self.beam_energy, 
            self.fragment.mass, 
            self.residual.mass, 
            fragment_angle * np.pi / 180
        )

        s = Reaction.__s_factor(
            self.beam.mass, 
            self.beam_energy, 
            self.fragment.mass, 
            self.residual.mass, 
            self.reaction_quit(residual_state)
        )

        return (r + np.sqrt(r ** 2 + s)) ** 2
    
    def residual_energy(self, residual_state: float) -> np.ndarray:
        r = Reaction.__r_factor(
            self.beam.mass,
            self.beam_energy, 
            self.residual.mass,
            self.fragment.mass,
            self.residual_angle(residual_state)
        )

        s = Reaction.__s_factor(
            self.beam.mass,
            self.beam_energy,
            self.residual.mass,
            self.fragment.mass,
            self.reaction_quit(residual_state)
        )

        return (r + np.sqrt(r ** 2 + s)) ** 2
    
    def residual_angle(self, residual_state: float, fragment_angle: float) -> float:
        fragment_ears = self.fragment_energy(residual_state)
        energy_relation = np.sqrt(self.beam.mass() * self.beam_energy / (self.fragment.mass() * fragment_ears))

        return np.pi / 2 - np.arctan(
            (energy_relation - np.cos(fragment_angle * np.pi / 180)) / np.sin(fragment_angle * np.pi / 180)
        )
    
    @staticmethod
    def __r_factor(beam_mass: float, beam_energy: float, 
                   instance_mass: float, partner_mass: float, angle: float) -> float:
        numerator = np.sqrt(beam_mass * instance_mass * beam_energy) * np.cos(angle)
        return numerator / (instance_mass + partner_mass)

    @staticmethod
    def __s_factor(beam_mass: float, beam_energy: float, 
                   instance_mass: float, partner_mass: float, reaction_quit: float) -> float:
        numerator = beam_energy * (partner_mass - beam_mass) + partner_mass * reaction_quit
        return numerator / (instance_mass + partner_mass)


if __name__ == '__main__':
    pass
