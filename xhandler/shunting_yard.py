from physics import *


class ReactionMaster:
    def __init__(self, notation: str) -> None:
        self.notation = notation

    def to_reaction(self) -> Reaction:
        pass

    def split_nucleus(self) -> list[Nuclei]:
        pass

    def to_nuclei(self, name: str) -> Nuclei:
        name = name.lower()
        
        nuclons = ReactionMaster.nuclons_from_name(name)
        charge = ReactionMaster.match_charge_by_name(name)
        return Nuclei(nuclons, charge)

    @staticmethod
    def nuclons_from_name(name: str) -> int:
        if name in ['p', 'd', 't']:
            return ['p', 'd', 't'].index(name) + 1

        pretend = ''
        for i in name:
            if i.isdigit():
                pretend += i

        return int(pretend)

    @staticmethod
    def match_charge_by_name(name: str) -> int:
        name = ''.join([char for char in name if char.isalpha()])
        match name:
            case 'h' | 'p' | 'd' | 't': return 1
            case 'he': return 2
            case 'li': return 3
            case 'be': return 4
            case 'b': return 5
            case 'c': return 6
            case 'n': return 7
            case 'o': return 8
            case 'f': return 9
            case _: return -1


if __name__ == '__main__':
    pass
