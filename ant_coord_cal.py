#!/bin/env python3
#
# Coordinate Transformations and Pointing Calibration
#
#   From the paper "Coordinate Transformations" Dr. John Stryjewski
#
#------------------------------------------------------------------------------

from skyfield.api import load, Star, wgs84, utc, Distance, Angle
import numpy as np
import datetime as D
import tools.find_peaks as find_peaks
from tools.antenna import Antenna
from scipy.optimize import least_squares as leastsq
from scipy.optimize import minimize

class Station():
    '''Info describing the ground station.
        antenna is an initiallized Antenna object
        decimal degrees and meters.
    '''
    #TODO: this class is rather redundant and could likely be fixed up
    def __init__(self, antenna, results):
        eph = load('de421.bsp')
        self.name = antenna.name
        self.lat = antenna.lat
        self.lon = antenna.lon
        self.alt = antenna.alt
        self.earth = eph['earth']
        self.antenna = wgs84.latlon(latitude_degrees=self.lat,
                                    longitude_degrees=self.lon,
                                    elevation_m=self.alt)
        self.ts = load.timescale()
        self.build_data_sets(results)

        print(self.antenna.latitude.dms(), self.antenna.longitude.dms(), self.antenna.elevation.m)

    def build_data_sets(self, results):
        ''' Convert our data sets into arrays of floats, in order
        '''
        # the radio sources canonical position as an ENU vector
        self.true = []
        # Where we observed the source as an ENU vector
        self.measured = []

        for row in results:
            # fill in our true and measured arrays.
            true_vect = self.true_enu_to_source(row[2])
            meas_vect = self.measured_enu_to_source((float(row[3]),float(row[4])))
           # print("\n--- Distance between in radians")
           # print(ENU().angle_between(true_vect, meas_vect))
            self.true.append(true_vect)
            self.measured.append(meas_vect)
        self.true = np.array(self.true)
        self.measured = np.array(self.measured)

    def measured_enu_to_source(self, angles):
        ''' ENU vector we measured to the source
        '''
        output = ENU(azel=angles)
       # print("\n--- measured angles")
       # print("{:.3f}, {:.3f}".format(output.az, output.el))
        return output.vector

    def true_enu_to_source(self, obs_time):
        ''' find the ENU vector to the source given an observation time. Hard
             Coded fro Cas_a for now.
        '''
        t = self.ts.utc(obs_time)
        self.src = Star(ra_hours=(23,23,24.000),
                        dec_degrees=(58,48,54.00))
        astrometric = (self.earth + self.antenna).at(t).observe(self.src).apparent()
        alt, az, rng = astrometric.altaz()
        output = ENU(azel=(az.degrees, alt.degrees))
       # print("\n--- true angles --- {}".format(str(obs_time)))
       # print("{:.3f}, {:.3f}".format(output.az, output.el))
        return output.vector

    def fit_data(self):
        ''' Try least squares fitting our data sets
        '''
        self.result = test_least_squares(self.measured, self.true)



class ENU():
    ''' Azimuth-Elevation-Range to East-North-Up frame
    '''
    def __init__(self, azel=None):
        if azel == None:
            self.from_enu(*[0,1,0])
        else:
            self.from_azel(*azel)
        self.north_vect = np.array([0, 1, 0])

    def from_azel(self, az, el, rng=1.0):
        ''' Calculate the ENU vector from an Azimuth, Elevation and Range
        '''
        self.az = az
        self.el = el
        self.rng = rng
        self.vector = np.array([np.sin(np.radians(self.az))*np.cos(np.radians(self.el)),
                                np.cos(np.radians(self.az))*np.cos(np.radians(self.el)),
                                np.sin(np.radians(self.el))], dtype=np.double).dot(rng)
        self.east = self.vector[0]
        self.north = self.vector[1]
        self.up = self.vector[2]

    def from_enu(self, e, n, u):
        ''' Build the vector and az el from the provided values
        '''
        self.vector = self.unit_vector(np.array([e, n, u]))
        self.east = self.vector[0]
        self.north = self.vector[1]
        self.up = self.vector[2]
        self.az = (np.degrees(np.arctan2(self.east, self.north)) + 360) % 360
        self.el = np.degrees(np.arcsin(self.up))
        self.rng = np.sqrt(e**2+n**2+u**2)

    def __str__(self):
        return "{}, {}, {}".format(str(self.vector), str(self.az), str(self.el))

    def unit_vector(self, vector):
        ''' returns the unit vector of the vector
        '''
        return vector / np.linalg.norm(vector)

    def angle_between(self, v1, v2):
        ''' Return the angle in radians between vectors 'v1' and 'v2'
        '''
        v1_u = self.unit_vector(v1)
        v2_u = self.unit_vector(v2)
        return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))


def R_east(theta):
    ''' Rotate around the east axis
    '''
    rads = np.radians(theta, dtype=np.double)
    matrix = np.array([[1, 0, 0],
                       [0, np.cos(rads, dtype=np.double), -np.sin(rads, dtype=np.double)],
                       [0, np.sin(rads, dtype=np.double), np.cos(rads, dtype=np.double)]], dtype=np.double)
    return matrix

def R_north(theta):
    ''' Rotate around the north axis
    '''
    rads = np.radians(theta, dtype=np.double)
    matrix = np.array([[np.cos(rads, dtype=np.double), 0, np.sin(rads, dtype=np.double)],
                       [0, 1, 0],
                       [-np.sin(rads, dtype=np.double), 0, np.cos(rads, dtype=np.double)]], dtype=np.double)
    return matrix

def R_up(theta):
    ''' Rotate around the up axis
    '''
    rads = np.radians(theta, dtype=np.double)
    matrix = np.array([[np.cos(rads, dtype=np.double), -np.sin(rads, dtype=np.double), 0],
                       [np.sin(rads, dtype=np.double), np.cos(rads, dtype=np.double), 0],
                       [0, 0, 1]], dtype=np.double)
    return matrix

def R_roll(psi):
    ''' Rotation in the roll axis, rotate the whole mount
    '''
    rads = np.radians(psi)
    matrix = np.array([[np.cos(rads), 0, np.sin(rads)],
                       [0, 1, 0],
                       [-np.sin(rads), 0, np.cos(rads)]])
    return matrix

def R_pitch(phi):
    ''' rotation in the pitch axis, rotate the whole mount
    '''
    rads = np.radians(phi)
    matrix = np.array([[1, 0, 0],
                       [0, np.cos(rads), -np.sin(rads)],
                       [0, np.sin(rads), np.cos(rads)]])
    return matrix

def R_yaw(theta):
    ''' Rotation in the yaw axis, rotate the whole mount
    '''
    rads = np.radians(theta)
    matrix = np.array([[np.cos(rads), np.sin(rads), 0],
                       [-np.sin(rads), np.cos(rads),0],
                       [0, 0, 1]])
    return matrix

def pointing_residual(variables, measured, true, printer=False):
    ''' The callback function for our fitting
    '''
    #TODO: Each of these commented out chunks are for more or less error terms.
#    az, el, non, skew, yaw, pitch, roll = variables
#    names = ['az', 'el', 'non', 'skew', 'yaw', 'pitch', 'roll']
#    model = (R_up(-az) @ R_north(non) @ R_east(el) @ R_north(-non) @ R_up(-skew) @ R_yaw(-yaw) @ R_pitch(-pitch) @ R_roll(-roll) @ true.T).T

#    az, el, non, yaw, pitch, roll = variables
#    names = ['az', 'el', 'non', 'yaw', 'pitch', 'roll']
#    model = (R_up(-az) @ R_north(non) @ R_east(el) @ R_north(-non) @ R_yaw(-yaw) @ R_pitch(-pitch) @ R_roll(-roll) @ true.T).T

#    yaw, pitch, roll = variables
#    names = ['yaw', 'pitch', 'roll']
#    model = (R_yaw(-yaw) @ R_pitch(-pitch) @ R_roll(-roll) @ true.T).T

    az, el, yaw, pitch, roll = variables
    names = ['az', 'el', 'yaw', 'pitch', 'roll']
    model = (R_up(-az) @ R_east(el) @ R_yaw(-yaw) @ R_pitch(-pitch) @ R_roll(-roll) @ true.T).T

#    az, el = variables
#    names = ['az', 'el']
#    model = (R_up(-az) @ R_east(el) @ true.T).T

    if printer == True:
        output = ""
        for i, value in enumerate(variables):
            output += "\n{}\t{:4f}".format(names[i], value)
        print(str(output))
        print("\n")
        rows, cols = model.shape
        count = 0
        while count < rows:
            enu_out = ENU()
            enu_out.from_enu(*measured[count,:])
            distance = np.degrees(enu_out.angle_between(model[count,:], measured[count,:]))
            output = ""
            output += "{:.4f}".format(distance)
            output += "\taz {:.4f}\tel {:.4f}".format(enu_out.az, enu_out.el)
            print(output)
            count += 1
        print("\n")
#    print(true.shape)

    return np.sum(np.subtract(1, np.sum(model*measured, axis=1)))

def test_least_squares(measured, true):
    ''' Try out some scipy least squares fitting
    '''
    # all of my "guesses" to start the fit
    azBias = 0.1
    elBias = 0.1
    az = -0.25
    non = 0.1
    el = 0.1
    skew = 0.1
    roll = 0.1
    pitch = 0.1
    yaw = 0.1
    droop = 0.1

    variables = [az, el, yaw, pitch, roll]
#    bounds = [az_bound, el_bound]
#    out = minimize(pointing_residual, variables, args=(measured, true),
#                   method='SLSQP',
#                   bounds=bounds,
#                   tol=1e-12)
    bounds = [-10.0, 10.0]
    out = leastsq(pointing_residual, variables, args=(measured, true),
                  bounds=bounds,
                  ftol=1e-12)
    return out


if __name__ == "__main__":
    tst_name = "real"
    if tst_name == "real":
        DATA_FILE_PATH = 'data/Band_Power1.csv'
        TIME_FILE_PATH = 'data/26M_CAS_A_samples.txt'
        data_set = find_peaks.get_csv_data(DATA_FILE_PATH)
        time_set = find_peaks.get_csv_data(TIME_FILE_PATH)
        results = find_peaks.test_times(time_set, 15, data_set)

        # for testing...
    #    print("\n{}\n".format(results[0]))
        antenna = Antenna(csvsource='tools/FCDAS_antennas.csv')
        antenna.ready('26M')
        station = Station(antenna, results)
        station.fit_data()
        print("\n")
        print(station.result)
        lst_out = pointing_residual(station.result.x, station.measured, station.true, printer=True)

    if tst_name == "test":
        #   test data
        measured = np.array([[0.17101007, 0.96984631, 0.17364818],
                             [0.17101007, 0.96984631, 0.17364818],
                             [0.17101007, 0.96984631, 0.17364818],
                             [0.17101007, 0.96984631, 0.17364818]], dtype=np.double)

        true = np.array([[0,1,0],
                         [0,1,0],
                         [0,1,0],
                         [0,1,0]], dtype=np.double)

        output = test_least_squares(measured, true)
        print(output)

        test_2 = (R_up(10.0) @ R_east(10.0) @ true.T).T
        pointing = ENU()
        pointing.from_enu(*test_2[1,:])
        print("\n--- measured offset for testing, this is the inital error.")
        print(pointing)

        print("\n--- lst out, this output should be 'corrected' back to 0")
        lst_out = pointing_residual(output.x, measured, true, printer=True)
        print(output.x)

        print("\n--- If lst was high precision, residual should be")
        lst_out = pointing_residual([-10.0, -10.0], measured, true, printer=True)
        print(lst_out)

        print("\n--- Rotate in and then out, should be close to 10s and 0s.")
        test_3 = (R_up(10.0) @ R_east(10.0) @ true.T).T
        pointing.from_enu(*test_3[1,:])
        print(pointing)
        # dont forget that order of operations matters, unwrap the rotations
        # by moving backwards through the steps.
        test_3 = (R_east(-10.0) @ R_up(-10.0) @ test_3.T).T
        pointing.from_enu(*test_3[1,:])
        print(pointing)

