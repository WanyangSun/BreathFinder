
import binascii
import struct
import time
import zlib
import xmltodict

# 定义谱图类
class Spectrum:
    def __init__(self, filename, scan, index, tic, peaks, mz, polarity, charge, ms_level, rt_in_second):
        self.filename = filename
        self.scan = scan
        self.tic = tic
        self.peaks = peaks
        self.mz = mz
        self.polarity = polarity
        self.charge = charge
        self.index = index
        self.ms_level = ms_level
        self.rt_in_second = rt_in_second
    def get_mgf_string(self):
        return "PLACE HOLDER"
    def get_mgf_peak_string(self):
        output_string = ""
        for peak in self.peaks:
            output_string += str(peak[0]) + "\t" + str(peak[1]) + "\n"
        return output_string
    @staticmethod
    def get_tsv_header():
        return "filename\tspectrumindex\tspectrumscan\tcharge\tmz"
    def get_max_mass(self):
        max_mass = 0.0
        for peak in self.peaks:
            max_mass = max(max_mass, peak[0])
        return max_mass
    #Straight up cosine between two spectra
#    def cosine_spectrum(self, other_spectrum, peak_tolerance):
#        total_score, reported_alignments = spectrum_alignment.score_alignment(self.peaks, other_spectrum.peaks, self.mz, other_spectrum.mz, peak_tolerance)
#        return total_score

    #def annotate_peaks(self, peak_annotations_list):
    #    for peak_annotation in peak_annotation_list:


#Decode peaks for mzXML
def decode_spectrum(line, peaks_precision, peaks_compression, struct_iter_ok):
    """https://groups.google.com/forum/#!topic/spctools-discuss/qK_QThoEzeQ"""
    decoded = binascii.a2b_base64(line)
    number_of_peaks = 0
    unpack_format1 = ""
    if peaks_compression == "zlib":
        decoded = zlib.decompress(decoded)
    #Assuming no compression
    if peaks_precision == 32:
        number_of_peaks = len(decoded)/4
        unpack_format1 = ">%df" % number_of_peaks
    else:
        number_of_peaks = len(decoded)/8
        unpack_format1 = ">%dd" % number_of_peaks
    # peaks = []
    # if struct_iter_ok:
    #     peak_iter = struct.iter_unpack(unpack_format1,decoded)
    #     peaks = [
    #        pair for pair in zip(*[peak_iter] * 2)
    #     ]
    # else:
    peaks = [
       pair for pair in zip(*[iter(struct.unpack(unpack_format1,decoded))] * 2)
    ]
    return peaks
    # peaks_list = struct.unpack(unpack_format1,decoded)
    # return [
    #     (peaks_list[i*2],peaks_list[i*2+1])
    #     for i in range(0,int(len(peaks_list)/2))
    # ]


def load_mzxml_file(filename):
    output_ms1 = []
    output_ms2 = []
    struct_iter_ok = True
    canary = True
    with open(filename) as fd:
        xmltodict_start = time.time()
        mzxml = xmltodict.parse(fd.read())
        xmltodict_end = time.time()
        print("mzXML file reading time: " + str(round((xmltodict_end - xmltodict_start), 4)) + ' s')
        read_scans = mzxml['mzXML']['msRun']['scan']
        filename_output = filename.split('\\')[-1]
        # filename_output = os.path.split(filename)[1]
        index = 1
        for index, scan in enumerate(read_scans):
            ms_level, spectrum, struct_iter_ok, canary = read_mzxml_scan(scan, index, filename_output, struct_iter_ok, canary)
            index += 1
            if ms_level == 1:
                output_ms1.append(spectrum)
            if ms_level == 2:
                output_ms2.append(spectrum)
            nested_scans = scan.get('scan',[])
            if not isinstance(nested_scans,list):
                nested_scans = [nested_scans]
            for nested_scan in nested_scans:
                ms_level, spectrum, struct_iter_ok, canary = read_mzxml_scan(nested_scan, index, filename_output, struct_iter_ok, canary)
                index += 1
                output_ms2.append(spectrum)
            # print(index)
    print('A total of %s MS1 spectra and %s MS2 spectra were extracted.' 
      % (len(output_ms1), len(output_ms2)))
    return output_ms1, output_ms2


def read_mzxml_scan(scan, index, filename_output, struct_iter_ok, canary):
    ms_level = int(scan['@msLevel'])
    scan_number = int(scan['@num'])
    rt_in_second = float(scan['@retentionTime'][2:-1])
    polarity = str(scan['@polarity'])
    tic = float(scan['@totIonCurrent'])
    #Optional fields
    base_peak_intensity = 0.0
    base_peak_mz = 0.0
    base_peak_intensity = float(scan.get('@basePeakIntensity', 0.0))
    base_peak_mz = float(scan.get('@basePeakMz', 0.0))
    try:
        precursor_mz_tag = scan['precursorMz']
        precursor_mz = float(precursor_mz_tag['#text'])
        precursor_scan = int(precursor_mz_tag.get('@precursorScanNum', 0))
        precursor_charge = int(precursor_mz_tag.get('@precursorCharge', 0))
        precursor_intensity = float(precursor_mz_tag.get('@precursorIntensity', 0))
    except:
        if ms_level == 2:
            raise
    peaks_precision = float(scan['peaks'].get('@precision', '32'))
    peaks_compression = scan['peaks'].get('@compressionType', 'none')
    peak_string = scan['peaks'].get('#text', '')
    if canary and peak_string != '':
        try:
            decode_spectrum(peak_string, peaks_precision, peaks_compression, struct_iter_ok)
        except:
            struct_iter_ok = False
        canary = False
    if peak_string != '':
        peaks = decode_spectrum(peak_string, peaks_precision, peaks_compression, struct_iter_ok)
    else:
        peaks = None
    if ms_level == 1:
        output = Spectrum(
            filename_output,
            scan_number,
            index,
            tic,
            peaks,
            0,
            polarity,
            0,
            ms_level,
            rt_in_second
        )
    if ms_level == 2:
        output = Spectrum(
            filename_output,
            scan_number,
            index,
            0,
            peaks,
            precursor_mz,
            polarity,
            precursor_charge,
            ms_level,
            rt_in_second
        )
    return ms_level, output, struct_iter_ok, canary
