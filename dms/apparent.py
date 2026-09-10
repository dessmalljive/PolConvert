import casacore.tables
from casacore.measures import measures
import astropy.time
from astropy import units as u 
import numpy as np

def get_apparent_source_positions(msfn):
    dm = measures()
    # First: get the reference time of the experiment:
    # (FIXME: I think we're supposed to round it to something)
    obst = casacore.tables.table(f'{msfn}::OBSERVATION')
    t0js = obst.getcol('TIME_RANGE')[0, 0]
    obst.close()
    t0jd = t0js/(u.d).to(u.s)
    t = astropy.time.Time(t0jd, format='mjd')
    expt_time = t.to_datetime() # as a datetime
    expt_epoch = dm.epoch('utc', str(expt_time))
    # Set our frame to that
    dm.do_frame(expt_epoch)
    field_t = casacore.tables.table(f"{msfn}::FIELD")
    # Following tConvert, we use the pointing centre as the direction to tweak:
    field_pos = field_t.getcol('DELAY_DIR')[:, 0, :]
    field_dict = field_t.getcolkeywords('DELAY_DIR')
    field_t.close()
    # Better to be intolerant than wrong, so we make demands on the units and epoch:
    assert field_dict['QuantumUnits'] == ['rad', 'rad']
    assert field_dict['MEASINFO']['Ref'] == 'J2000'
    app_poses = np.zeros_like(field_pos)
    for i in range(field_pos.shape[0]):
        ra, dec = field_pos[i] # in radians or else
        svec = dm.direction('j2000', f'{ra}rad', f'{dec}rad')
        # And find the adjusted position:
        app_pos = dm.measure(svec, 'app')
        assert app_pos['m0']['unit'] =='rad' and app_pos['m1']['unit'] =='rad'
        app_ra = app_pos['m0']['value']
        app_dec = app_pos['m1']['value']                
        app_poses[i, 0] = app_ra
        app_poses[i, 1] = app_dec
    return app_poses

if __name__=='__main__':
    msfn = 'em155a-scratch.ms'
    app = get_apparent_source_positions(msfn)
    print(app)


