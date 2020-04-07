#!/usr/bin/python3 -u

from collections import namedtuple


class RMSD_Hit(namedtuple('RMSD_Hit', ['pdb', 'chain', 'model', 'res_start', 'res_end', 'hit_sequence', 'rmsd'])):
    def __str__(self):
        out = "RMSD_HIT: {} model= {:>3}, chain= {:1} hit= {:>4} {} {:<4} rms= {:>6.4f}"
        out = out.format(self.pdb, self.model, self.chain,
                         self.res_start, self.hit_sequence, self.res_end, self.rmsd)
        return out


class Water_Hit(namedtuple('Water_Hit', ['water_id', 'water_rms'])):
    def __init__(self):
        super().__init__('Water_Hit')
        self.water_id = None
        self.water_rms = None

    def __str__(self):
        str = ""
        if self.water_id and self.water_rms:
            str = "WAT= {:4} wat_rms= {:>6.4f}"
            str.format(format(self.water_id, self.water_rms))
        return str
