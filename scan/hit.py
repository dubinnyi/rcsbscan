#!/usr/bin/python3 -u

from collections import namedtuple


class Hit(namedtuple('RMSD_Hit', ['pdb', 'chain', 'model', 'res_start', 'res_end', \
                                  'hit_sequence', 'rmsd'])):
    def add_Water(self, water_id, water_rms):
        self.water_id = water_id
        self.water_rms = water_rms

    def add_Structure(self, structure, method='', chain_size=None):
        self.structure = structure
        self.method = method
        self.chain_size=chain_size

    def __str__(self):
        out = "{} {} model= {:>3}, chain= {:1} size= {:>4} hit= {:>4} {} {:<4} rms= {:>6.4f}"
        out = out.format(self.pdb, self.method, self.model, self.chain, self.chain_size,
                         self.res_start, self.hit_sequence, self.res_end, self.rmsd)
        if hasattr(self, 'water_id'):
            out += " WAT= {:4} wat_rms= {:>6.4f}"
            out = out.format(self.water_id, self.water_rms)
        return out



#class Water_Hit(namedtuple('Water_Hit', ['water_id', 'water_rms'])):
#    def __str__(self):
#        out = ""
#        if self.water_id and self.water_rms:
#            out = "WAT= {:4} wat_rms= {:>6.4f}"
#            out = out.format(self.water_id, self.water_rms)
#        return out
# </editor-fold>

# class Hit(RMSD_Hit, Water_Hit):
#    def __init__
