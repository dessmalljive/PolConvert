import numpy as np

tb.open('em155a-scratch.ms')
t2  = tb.query('SCAN_NUMBER=1 && ANTENNA="T6"')

c = t2.getcol('CORRECTED_DATA')[0][0,0]
print(f"{c=}")
