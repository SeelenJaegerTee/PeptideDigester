import os
import datetime
from io import *
import logging


class Config:
    DEBUG = True
    CHUNK_SIZE = 65535      # 2^16 -1 just in case


class Setup:
    __already_executed = False

    def __init__(self, input_path='.', scan_subfolders=False):
        if Setup.__already_executed is False:
            Setup.files_index = []
            input_path.replace('\\', '/')
            if input_path == '.' or input_path == '':
                input_path = (os.getcwd()+'\\').replace('\\', '/')
            elif input_path[-1] != '/':
                input_path += '/'

            def scan_folder(path):              # inner function because ...
                for f in os.listdir(path):
                    if os.path.isfile(path+f):
                        Setup.files_index.append(path+f)
                    elif f != '.idea' and f != 'output':
                        if scan_subfolders is True:         # ... we may want to do recursion ...
                            scan_folder(path+f+'/')

            scan_folder(input_path)             # ... at this point
            logging.basicConfig(filename='logfile.log', level=logging.DEBUG)
            logging.info('Setup initiated: '+str(datetime.datetime.now()))
            Setup.__already_executed = True
        else:
            logging.warning('warning, tried to execute setup twice')

    def extract_filetype(self, file_ending):
        ret = []
        for f in Setup.files_index:
            if file_ending in f:
                ret.append(f)
        return ret


# can provide lines from the file by buffering a large chunk of the file
# should prevent filling memory in case input file is incredibly large
# warning when requesting a line and if the file has no line-breaks (\n)
# the entire file will be successively stored into a string
# at this point it would probably have been faster to just readlines the entire file
class ChunkProvider:
    def __init__(self, file_path):
        self.file_path = file_path
        self.pos = 0
        self.chunk_buffer = ''
        self.EOF_reached = False
        self.file_exhausted = False

        self.line = 0
        self.chunk_pos = 0
        self.read_chunk(Config.CHUNK_SIZE)

    def read_chunk(self, chunk_size=Config.CHUNK_SIZE):
        if Config.DEBUG is True:
            logging.debug('reading chunk from position '+str(self.pos)+' onward ('+str(Config.CHUNK_SIZE)+') bytes')
        if not self.EOF_reached:
            with open(self.file_path) as fp:
                fp.seek(self.pos)
                self.chunk_buffer = fp.read(chunk_size)
                # if Config.DEBUG is True: print('chunk loaded: ' + repr(self.chunk_buffer))
                self.pos = fp.tell()
                # check whether EOF is reached
                if self.pos == fp.seek(0, SEEK_END):
                    self.EOF_reached = True
                    if Config.DEBUG is True:
                        logging.debug('EOF reached')
                # we need to reset the filepointer to where it was
                fp.seek(self.pos)
        else:
            self.chunk_buffer = ''
            if Config.DEBUG is True:
                logging.debug('can not load further chunks, reached end of file. Returning empty string')

    def next_line(self):
        if Config.DEBUG:
            logging.debug('trying to retrieve next line')
        new_line_pos = self.chunk_buffer.find('\n', self.chunk_pos)
        if new_line_pos != -1:
            ret = self.chunk_buffer[self.chunk_pos: new_line_pos+1]
            self.line += 1
            if Config.DEBUG:
                logging.debug('found line at pos: '+str(self.pos))
            self.chunk_pos = new_line_pos + 1
        elif not self.EOF_reached:
            if Config.DEBUG:
                logging.debug('found no newline but there is still more file to buffer ...')
            ret = self.chunk_buffer[self.chunk_pos:]
            self.read_chunk()
            self.chunk_pos = 0
            ret += self.next_line()
        elif not self.file_exhausted:
            if Config.DEBUG:
                logging.debug('last bit of file already loaded in buffer retrieving until EOF')
            ret = self.chunk_buffer[self.chunk_pos:]
            self.line += 1
            self.file_exhausted = True
        else:
            if Config.DEBUG:
                logging.debug('finished, file contents are exhausted. Returning empty string')
            return ''
        return ret


class File:
    def __init__(self, file_path, name):
        self.name = name
        self.cp = ChunkProvider(file_path)
        self.peptides = []

    def analyze(self):
        buff = self.cp.next_line()
        while buff != '':
            if buff[0] is '>':
                self.peptides.append(Peptide(buff.rstrip()))
            else:
                self.peptides[-1].sequence = buff
            buff = self.cp.next_line()
        for pep in self.peptides:
            pep.analyze()


class Peptide:
    def __init__(self, name):
        self.name = name.rstrip()
        self.amino_acids = {}
        self.sequence = ''
        self.restriction_sites = []

    def analyze(self):
        for aa in self.sequence:
            if aa in self.amino_acids:
                self.amino_acids[aa] += 1
            else:
                self.amino_acids[aa] = 1


class Protease:
    def __init__(self, name, target_sequence, cleavage_pos):
        self.name = name
        self.target_sequence = target_sequence
        self.cleavage_pos = cleavage_pos


class RestrictionSite:
    def __init__(self, position, protease):
        self.position = position
        self.protease = protease


#########################################
# Script
#########################################
files_index = Setup(scan_subfolders=False).extract_filetype('.fasta')
print(files_index)
files = []
for f in files_index:
    buff = ''
    for c in f[::-1]:
        buff += c
        if c == '/':
            buff = buff[:-1]
            break
    buff = buff[::-1]
    print(buff)

    files.append(File(f, buff))

print(files)
for f in files:
    print(f.name)
    f.analyze()
    for pep in f.peptides:
        print(pep.name)
        print(pep.sequence)
