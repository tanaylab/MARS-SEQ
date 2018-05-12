#!/usr/bin/env python2
import sys
import collections
import yaml
import argparse
import pandas as pd
import barcodes_v1
import barcodes_v3


########################################################################
ROWS = [chr(c+ord('A')) for c in xrange(16)]
COLS = [str(c+1) for c in xrange(24)]


########################################################################
Plate = collections.namedtuple('Plate', ['user',
                                         'plate_id',
                                         'subject_id',
                                         'amp_batch',
                                         'barcodes',
                                         'spike_type',
                                         'spike_dilution',
                                         'spike_volume',
                                         'empty_wells'])


########################################################################
def main(argv):
    args = parse_args(argv[1:])

    plates = read_plates(args.plates)

    wells = []
    for plate in plates:
        wells.append(plate_to_wells(plate))
    wells = pd.concat(wells, ignore_index=True)

    wells.to_csv(args.output, sep='\t', header=True, index=False)


########################################################################
def parse_args(argv):
    parser = argparse.ArgumentParser(description = "Generate the mars-seq wells_cells annotation.")
    parser.add_argument('plates', metavar='plates.yaml', type=str,
                        help="YAML file holding plates configuraion.")
    parser.add_argument('output', metavar='output.txt', type=str, nargs='?', default='wells_cells.txt',
                        help="Output filename")

    return parser.parse_args(argv)


########################################################################
def read_plates(plates_fn):
    config = yaml.load(open(plates_fn, 'r'))

    defaults = config['defaults']
    plates = config['plates']

    for id in plates:
        plates[id]['plate_id'] = str(id)

    plates = plates.values()
    def_keys = frozenset(defaults.keys())
    for plate in plates:
        missing = def_keys - frozenset(plate.keys())
        for key in missing:
            plate[key] = defaults[key]
        plate['empty_wells'] = plate['empty_wells'].split()
        if (isinstance(plate['amp_batch'], str)):
            plate['amp_batch'] = [plate['amp_batch']]
        if (len(plate['amp_batch']) == 1):
            plate['amp_batch'] *= 2

    plates = [Plate(**p) for p in plates]

    plates.sort(key=lambda x: x.plate_id)

    return plates


########################################################################
def plate_to_wells(plate):
    barcodes = choose_barcodes(plate.barcodes)
    if (barcodes is None):
        raise ValueError("Invalid barcode identifier '%s' in plate '%s'." % (plate.barcodes, plate.plate_id))
    wells = []
    for col_shift in xrange(2):
        for row_shift in xrange(2):
            qp = quarter_plate_to_wells(plate, row_shift, col_shift, barcodes[row_shift][col_shift])
            wells.append(qp)
    wells = pd.concat(wells, ignore_index=True)
    wells.Well_ID = ['W%s%s%03d' % (plate.user, plate.plate_id, c+1) for c in range(len(wells))]
    return wells


########################################################################
def choose_barcodes(barcodes_id):
    if (barcodes_id == 'v1_g1g2'):
        return ( (barcodes_v1.gr1, barcodes_v1.gr1),
                 (barcodes_v1.gr2, barcodes_v1.gr2) )
    if (barcodes_id == 'v3_g1g2'):
        return ( (barcodes_v3.gr1, barcodes_v3.gr1),
                 (barcodes_v3.gr2, barcodes_v3.gr2) )
    if (barcodes_id == 'v3_g1g3'):
        return ( (barcodes_v3.gr1, barcodes_v3.gr1),
                 (barcodes_v3.gr3, barcodes_v3.gr3) )
    if (barcodes_id == 'v3_g1g4'):
        return ( (barcodes_v3.gr1, barcodes_v3.gr1),
                 (barcodes_v3.gr4, barcodes_v3.gr4) )
    if (barcodes_id == 'v3_g2g3'):
        return ( (barcodes_v3.gr2, barcodes_v3.gr2),
                 (barcodes_v3.gr3, barcodes_v3.gr3) )
    if (barcodes_id == 'v3_g2g4'):
        return ( (barcodes_v3.gr2, barcodes_v3.gr2),
                 (barcodes_v3.gr4, barcodes_v3.gr4) )
    if (barcodes_id == 'v3_g3g4'):
        return ( (barcodes_v3.gr3, barcodes_v3.gr3),
                 (barcodes_v3.gr4, barcodes_v3.gr4) )
    if (barcodes_id == 'v3'):
        return ( (barcodes_v3.gr1, barcodes_v3.gr2),
                 (barcodes_v3.gr3, barcodes_v3.gr4) )
    return None

########################################################################
def quarter_plate_to_wells(plate, row_shift, col_shift, barcodes):
    coordinates = []
    cell_barcode = []
    number_of_cells = []

    for col in xrange(12):
        for row in xrange(8):
            row_str = ROWS[row*2+row_shift]
            col_str = COLS[col*2+col_shift]
            coord = row_str + col_str
            alt_coord = col_str + row_str
            empty = (coord in plate.empty_wells) or (alt_coord in plate.empty_wells)
            coordinates.append(coord)
            cell_barcode.append(barcodes[row][col])
            number_of_cells.append(0 if empty else 1)


    df = collections.OrderedDict()
    df['Well_ID'] = 'xxx'
    df['well_coordinates'] = coordinates
    df['plate_ID'] = 'P%s%s' % (plate.user, plate.plate_id)
    df['Subject_ID'] = plate.subject_id
    df['Amp_batch_ID'] = plate.amp_batch[row_shift]
    df['Cell_barcode'] = cell_barcode
    df['Spike_type'] = plate.spike_type
    df['Spike_dilution'] = ("%.10f" % plate.spike_dilution).rstrip('0')
    df['Spike_volume_ul'] = ("%.10f" % plate.spike_volume).rstrip('0')
    df['Number_of_cells'] = number_of_cells
    df['is_primer_added'] = 1
    df['Raw_RNA_amount'] = 'NA'
    df['Lysis_date'] = 'NA'
    df = pd.DataFrame(df)

    return df


########################################################################
def write_output(output_fn, plates):
    headers = [ 'Well_ID',
                'well_coordinates',
                'plate_ID',
                'Subject_ID',
                'Amp_batch_ID',
                'Cell_barcode',
                'Spike_type',
                'Spike_dilution',
                'Spike_volume_ul',
                'Number_of_cells',
                'is_primer_added',
                'Raw_RNA_amount',
                'Lysis_date' ]
    _rows = 'ABCDEFGHIJKLMNOP'
    _cols = [str(c) for c in xrange(1,25)]


########################################################################
if (__name__ == '__main__'):
    sys.exit(main(sys.argv))


# class _PlateParams:
#     headers = [ 'Well_ID',
#                 'well_coordinates',
#                 'plate_ID',
#                 'Subject_ID',
#                 'Amp_batch_ID',
#                 'Cell_barcode',
#                 'Spike_type',
#                 'Spike_dilution',
#                 'Spike_volume_ul',
#                 'Number_of_cells',
#                 'is_primer_added',
#                 'Raw_RNA_amount',
#                 'Lysis_date' ]
#     _rows = 'ABCDEFGHIJKLMNOP'
#     _cols = [str(c) for c in xrange(1,25)]
#
#     def __init__(self, params):
#         self.subject        = self._getStrParam(params, 'subject')
#         self.amp_batch_a    = self._getStrParam(params, 'amp_batch_a')
#         self.amp_batch_h    = self._getStrParam(params, 'amp_batch_h')
#         self.barcodes       = self._getStrParam(params, 'barcodes')
#         self.spike_type     = self._getStrParam(params, 'spike_type')
#         self.spike_dilution = self._getStrParam(params, 'spike_dilution')
#         self.spike_volume   = self._getStrParam(params, 'spike_volume')
#         self.empty_wells    = self._getStrParam(params, 'empty_wells')
#
#         user_initials = self._getStrParam(params, 'user_initials')
#         plate_no      = self._getIntParam(params, 'plate_no')
#         if ((user_initials == "#ERROR") or (plate_no == "#ERROR")):
#             self.well_prefix = '#ERROR'
#             self.plate       = '#ERROR'
#         else:
#             self.well_prefix = 'W%s%04d' % (user_initials, plate_no)
#             self.plate       = 'P%s%04d' % (user_initials, plate_no)
#
#         self.empty_wells = self.empty_wells.split()
#         self.empty_wells = set([w if w[0].isalpha() else w[-1]+w[:-1] for w in self.empty_wells])
#
#     def cellPos(self, row, col):
#         return self._rows[row] + self._cols[col]
#
#     def _getStrParam(self, params, name):
#         param = params.get(name)
#         if ((param is None) or (param == '')):
#             return "#ERROR"
#         return param
#
#     def _getIntParam(self, params, name):
#         param = self._getStrParam(params, name)
#         if ((param is None) or (not param.isdigit())):
#             return "#ERROR"
#         return int(param)
#
#     def _getBoolParam(self, params, name):
#         param = self._getStrParam(params, name)
#         if (param == 'on'):
#             return True
#         return False
#
#
# class WellsCells(webapp2.RequestHandler):
#
#     def post(self):
#         self.response.headers['Content-Type'] = 'text/plain'
#         self.response.headers['Content-Disposition'] = 'attachment; filename="wells_cells.txt"'
#         self.response.write('\t'.join(_PlateParams.headers) + '\n')
#         # Extract post items
#         items = self.request.POST.items()
#         items = [(key.encode('ascii', 'ignore'), value.encode('ascii', 'ignore')) for (key, value) in items]
#         # Get globals (all items with a '_0' suffix)
#         suffix_pos = -2
#         suffix = '_0'
#         globs = {key[:suffix_pos]:value for (key, value) in items if (key[suffix_pos:] == suffix)}
#         items = [(key, value) for (key, value) in items if (key[suffix_pos:] != suffix)]
#         # Go over plate parameters and write plates
#         while (len(items) > 0):
#             id = items[0][0]
#             suffix_pos = id.rindex('_')
#             if (suffix_pos < 0):
#                 break
#             suffix_pos -= len(id)
#             suffix = id[suffix_pos:]
#             fields = {key[:suffix_pos]:value for (key, value) in items if (key[suffix_pos:] == suffix)}
#             items = [(key, value) for (key, value) in items if (key[suffix_pos:] != suffix)]
#             fields.update(globs)
#             self.printPlate(_PlateParams(fields))
#
#     def printPlate(self, plate_params):
#         barcodes = self.choose_barcodes(plate_params)
#         if (plate_params.amp_batch_a != "#ERROR"):
#             self.printQuarterPlate(plate_params, 0, 0, plate_params.amp_batch_a, barcodes[0][0], 1)
#             self.printQuarterPlate(plate_params, 1, 0, plate_params.amp_batch_a, barcodes[1][0], 1 + 96)
#         if (plate_params.amp_batch_h != "#ERROR"):
#             self.printQuarterPlate(plate_params, 0, 1, plate_params.amp_batch_h, barcodes[0][1], 1 + 96*2)
#             self.printQuarterPlate(plate_params, 1, 1, plate_params.amp_batch_h, barcodes[1][1], 1 + 96*3)
#
#
#     def choose_barcodes(self, plate_params):
#         if (plate_params.barcodes == 'v1 GR1-GR2'):
#             return ( (barcodes_v1.gr1, barcodes_v1.gr1),
#                      (barcodes_v1.gr2, barcodes_v1.gr2) )
#         if (plate_params.barcodes == 'v3 GR1-GR2'):
#             return ( (barcodes_v3.gr1, barcodes_v3.gr1),
#                      (barcodes_v3.gr2, barcodes_v3.gr2) )
#         if (plate_params.barcodes == 'v3 GR1-GR3'):
#             return ( (barcodes_v3.gr1, barcodes_v3.gr1),
#                      (barcodes_v3.gr3, barcodes_v3.gr3) )
#         if (plate_params.barcodes == 'v3 GR1-GR4'):
#             return ( (barcodes_v3.gr1, barcodes_v3.gr1),
#                      (barcodes_v3.gr4, barcodes_v3.gr4) )
#         if (plate_params.barcodes == 'v3 GR2-GR3'):
#             return ( (barcodes_v3.gr2, barcodes_v3.gr2),
#                      (barcodes_v3.gr3, barcodes_v3.gr3) )
#         if (plate_params.barcodes == 'v3 GR2-GR4'):
#             return ( (barcodes_v3.gr2, barcodes_v3.gr2),
#                      (barcodes_v3.gr4, barcodes_v3.gr4) )
#         if (plate_params.barcodes == 'v3 GR3-GR4'):
#             return ( (barcodes_v3.gr3, barcodes_v3.gr3),
#                      (barcodes_v3.gr4, barcodes_v3.gr4) )
#         if (plate_params.barcodes == 'v3 GR1-2-3-4'):
#             return ( (barcodes_v3.gr1, barcodes_v3.gr2),
#                      (barcodes_v3.gr3, barcodes_v3.gr4) )
#         return None
#
#
# app = webapp2.WSGIApplication([
#     ('/wells_cells.txt', WellsCells),
# ], debug=True)
