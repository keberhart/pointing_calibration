#!/usr/bin/env python
#
# create a list of points so we can track a target with the new Orbit ACU
#
#   9APR19 - Kyle Eberhart
#   5MAY20 - Kyle Eberhart - extend this to other targets and antennas
#
#   This project made use of Skyfield, http://rhodesmill.org/skyfield/
#
#----------------------------------------------------------------------

from skyfield.api import Star, load, wgs84, position_of_radec
from datetime import timedelta
from numpy import arange

class Target():

    def __init__(self):
        # load the planetary positions and our position
        self.planets = load('de421.bsp')
        self.earth = self.planets['earth']
        self.Cas_A = Star(ra_hours=(23, 23, 24.000),
                          dec_degrees=(58, 48, 54.00))
        # The coords below were pulled from Guy Kauffmans notes on the
        #  Old 26M OCU. They don't track like I would expect...
#        self.Cas_A = Star(ra_hours=(23, 23, 5.00),
#                           dec_degrees=(58, 46, 0.00),
#                           epoch=1992.5)

        self.sources = {"CAS_A":self.Cas_A,
                        "SUN":self.planets['sun'],
                        "MOON":self.planets['moon']}

    def get_list(self):
        return list(self.sources)

    def ready_antenna(self, lat, lon, alt):
        self.ant = self.earth + wgs84.latlon(latitude_degrees=lat,
                        longitude_degrees=lon,
                        elevation_m=alt)

    def ready(self, target):
        self.name = self.sources[target]

    def find_angles(self, cust_time, peak_az, peak_alt):
        ''' find angles to our source. use a custom time, so we can find
            differences in pointing.
            provide a tuple of (year, month, day, hours, min, sec)
        '''
        ts = load.timescale(builtin=True)
        cust_time = ts.utc(*cust_time)
        astrometric = self.ant.at(cust_time).observe(self.name).apparent()
        alt, az, distance = astrometric.altaz()
        delta_az = az.degrees - peak_az
        delta_alt = alt.degrees - peak_alt
        output = []
        output.append('az_peak\t\taz_pred\t\tdelta\n')
        output.append('\t'.join([str(peak_az), str(az.degrees), str(delta_az)]))
        output.append('\n')
        output.append('\t'.join([str(peak_alt), str(alt.degrees), str(delta_alt)]))
        output.append('\n')
        print(self.ant)
        print(astrometric.position)
        special = self.ant.at(cust_time).from_altaz(alt_degrees=peak_alt, az_degrees=peak_az)
        print(special.position)
        print(astrometric.separation_from(special))


        for item in output:
            print(item)

    def drift_scan(self):
        ''' Create a list of points for drift scanning across the source.
            separate the points by DWELL_TIME minutes and DELTA degree in declination
            angle.
        '''
        DWELL_TIME = 15 # in minutes
        DELTA = 0 # declination offset in degrees
        DURATION = 300 # in minutes
        ts = load.timescale(builtin=True)
        self.drift_script = []
        self.drift_samples = []
        split_int = DWELL_TIME
        split_time = timedelta(minutes=split_int)
        sample_time = ts.now()+timedelta(minutes=15)
        move_time = sample_time - timedelta(minutes=(split_int/2))
        # define the offset, but maybe make it a variable later
        offset = DELTA
        # figure out the maximum declination angle
        max_dec = self.name.dec.degrees + (offset * 3)
        # set our starting point
        self.name.dec = self.name.dec.degrees - (offset * 3)
        # set up another breakout method for duration of test.
        max_samples = int(DURATION / DWELL_TIME)
        curr_sample = 0
        while self.name.dec <= max_dec:
            if curr_sample > max_samples:
                break
            # find angles at time
            astrometric = self.ant.at(sample_time).observe(self.name).apparent()
            alt, az, distance = astrometric.altaz()
            t_str = sample_time.utc_strftime('%Y:%m:%d:%H:%M:%S')
            m_str = move_time.utc_strftime('%Y-%j-%H:%M:%S')
            az_str = "{:.4f}".format(az.degrees)
            alt_str = "{:.4f}".format(alt.degrees)

            self.drift_script.append('wait @(time("{}"))\n'.format(m_str))
            self.drift_script.append('start point_26M({}, {})\n'.format(az_str, alt_str))
            self.drift_script.append('write "Peak @ {}"\n'.format(t_str))
            self.drift_samples.append("{},{},{}\n".format(t_str, az_str, alt_str))
#            self.drift_list.append(', '.join([m_str, t_str, az_str, alt_str,'\n']))
            # increment
            sample_time = sample_time + split_time
            move_time = move_time + split_time
            self.name.dec = self.name.dec + offset
            curr_sample = curr_sample + 1

    def generate_report(self):
        # the tricky bits to get look angles for a list of times
        ts = load.timescale(builtin=True)

        # get rid of the milliseconds in there
        truth = ts.now().utc
        now_ish = ts.utc(truth[0],truth[1],truth[2],truth[3],truth[4])
        now = now_ish.tt

        # make an arry with times incremented by really close to 30 sec
        # this looks weird since it is terrestrial time floating point
        run_time = arange(now, now + .5, .000347222)

        # feed that array to the time function as julian days
        t = ts.tt(jd=run_time)

        # spit out the positon and vector of our location
        astrometric = self.ant.at(t).observe(self.name).apparent()

        # convert this to alt and az lists
        alt, az, distance = astrometric.altaz()

        t_list = t.utc_strftime('%Y:%m:%d:%H:%M:%S')
        alt_list = alt.degrees.tolist()
        az_list = az.degrees.tolist()

        self.report = []
        self.rise_time = None
        self.fade_time = None

        for idx,item in enumerate(t_list):
            if alt_list[idx] < 10.0:
                continue
            else:
                if self.rise_time is None:
                    self.rise_time = t_list[idx]
                self.fade_time = t_list[idx]
                az_out = "{:.4f}".format(az_list[idx])
                alt_out = "{:.4f}".format(alt_list[idx])
                output = t_list[idx]+", "+az_out+", "+alt_out+"\n"
                self.report.append(output)


