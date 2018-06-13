import os
import datetime
from io import *


class Config:
    DEBUG = True
    CHUNK_SIZE = 65535


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
        self.read_chunk()

    def read_chunk(self, chunk_size=Config.CHUNK_SIZE):
        # if Config.DEBUG: print('before call self.chunk_buffer', repr(self.chunk_buffer))
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
                        print('EOF reached')
                # we need to reset the filepointer to where it was
                fp.seek(self.pos)
        else:
            self.chunk_buffer = ''
            if Config.DEBUG is True:
                print('requested to read chunk despite EOF reached, returning empty string')
        # if Config.DEBUG: print('after call self.chunk_buffer', repr(self.chunk_buffer))

    def next_line(self):
        new_line_pos = self.chunk_buffer.find('\n', self.chunk_pos)
        if Config.DEBUG: print('new_line_pos', new_line_pos)
        if new_line_pos != -1:
            ret = self.chunk_buffer[self.chunk_pos: new_line_pos+1]
            self.line += 1
            if Config.DEBUG: print('found line')
            self.chunk_pos = new_line_pos + 1
        elif not self.EOF_reached:
            if Config.DEBUG: print('still more file to buffer ...')
            ret = self.chunk_buffer[self.chunk_pos: ]
            if Config.DEBUG: print('reloading')
            self.read_chunk()
            self.chunk_pos = 0
            ret += self.next_line()
        elif not self.file_exhausted:
            if Config.DEBUG: print('last bit of file already loaded in buffer')
            ret = self.chunk_buffer[self.chunk_pos:]
            self.line += 1
            self.file_exhausted = True
        else:
            if Config.DEBUG: print('finished')
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


my_file = File('C:/Users/NeugebauerC/Desktop/python kurs/Digester/ins_human.fasta', 'ins_human.fasta')
my_file.analyze()
print(my_file.name)
print(my_file.peptides[0].name)
print(my_file.peptides[0].sequence)
print(my_file.peptides[0].amino_acids)