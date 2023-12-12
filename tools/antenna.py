#/usr/bin/env python3

import csv

class Antenna():
    '''Interact with an antenna CSV file.'''

    def __init__(self, csvsource='FCDAS_antennas.csv'):
        '''open the file and put it in a dictionary.'''
        self.ants = {}
        with open(csvsource, newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                self.ants.update({row['Antenna']:row})

    def ready(self, name):
        '''read out the values for the selected antenna'''
        self.name = name
        self.source = self.ants[name]
        self.lat = float(self.source['DecLat'])
        self.lon = float(self.source['DecLon'])
        self.alt = float(self.source['DecAlt'])

    def get_list(self):
        return list(self.ants)

    def __str__(self):
        pile = []
        for row in self.ants:
            pile.append(row)
        return str(pile)


if __name__=='__main__':
    testing = Antenna()
    print(testing)
    print(testing.key_list())

