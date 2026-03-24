import requests
import re

# API URL for Horizons

url = 'https://ssd.jpl.nasa.gov/api/horizons.api'

def get_reference_ephemerides():

# So far I'm going to use Jupiter's actual center not it's system
# Barycenter, but I may have to change that later if it turns out to be 
# too inaccurate.
# The ranges for the time are tentative.

    #Io query

    params_io = {
        "format": "json",
        "COMMAND": "'501'", # Io ID
        "MAKE_EPHEM": "'YES'",
        "OBJ_DATA": "'NO'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@599'", # Jupiter Center
        "START_TIME": "'1600-01-11 00:00:00'",
        "STOP_TIME": "'2200-01-09 00:00:00'",
        "STEP_SIZE": "'10 d'",
        "VEC_TABLE": "'2'",
        "VEC_CORR": "'NONE'",
        "OUT_UNITS": "'KM-S'",
        "TIME_TYPE": "'TDB'",
        "CAL_FORMAT": "'JD'",
        "CSV_FORMAT": "'NO'",
        "VEC_LABELS": "'YES'",
    }

    r_io = requests.get(url, params=params_io, timeout=30)
    r_io.raise_for_status()
    data_io = r_io.json()

    if "error" in data_io:
        raise RuntimeError(data_io["error"])
    
    text_io = data_io["result"]

    block_io = text_io.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()

    num = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:E[+-]?\d+)?" # Maybe I should define this earlier

    records_io = re.split(r"(?=\s*\d+\.\d+\s*=)", block_io)

    vectors_io = []
    for rec in records_io:
        m = re.search(
            rf"(?P<jd>{num}).*?"
            rf"X\s*=\s*(?P<x>{num}).*?"
            rf"Y\s*=\s*(?P<y>{num}).*?"
            rf"Z\s*=\s*(?P<z>{num}).*?"
            rf"VX\s*=\s*(?P<vx>{num}).*?"
            rf"VY\s*=\s*(?P<vy>{num}).*?"
            rf"VZ\s*=\s*(?P<vz>{num})",
            rec,
            re.DOTALL,
        )
        if m:
            vectors_io.append([
                float(m.group("x")),
                float(m.group("y")),
                float(m.group("z")),
                float(m.group("vx")),
                float(m.group("vy")),
                float(m.group("vz")),
            ])

    print(vectors_io[0])
    print(len(vectors_io))

# This is just the step size for entries in the queried
# ephemeris. Our numerical methods will have to use much 
# smaller step sizes.

# Europa query

    params_europa = {
        "format": "json",
        "COMMAND": "'502'", # Europa ID
        "MAKE_EPHEM": "'YES'",
        "OBJ_DATA": "'NO'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@599'", # Jupiter Center
        "START_TIME": "'1600-01-11 00:00:00'",
        "STOP_TIME": "'2200-01-09 00:00:00'",
        "STEP_SIZE": "'10 d'",
        "VEC_TABLE": "'2'",
        "VEC_CORR": "'NONE'",
        "OUT_UNITS": "'KM-S'",
        "TIME_TYPE": "'TDB'",
        "CAL_FORMAT": "'JD'",
        "CSV_FORMAT": "'NO'",
        "VEC_LABELS": "'YES'",
    }

    r_europa = requests.get(url, params=params_europa, timeout=30)
    r_europa.raise_for_status()
    data_europa = r_europa.json()

    if "error" in data_europa:
        raise RuntimeError(data_europa["error"])

    text_europa = data_europa["result"]

    block_europa = text_europa.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()

    records_europa = re.split(r"(?=\s*\d+\.\d+\s*=)", block_europa)

    vectors_europa = []
    for rec in records_europa:
        m = re.search(
            rf"(?P<jd>{num}).*?"
            rf"X\s*=\s*(?P<x>{num}).*?"
            rf"Y\s*=\s*(?P<y>{num}).*?"
            rf"Z\s*=\s*(?P<z>{num}).*?"
            rf"VX\s*=\s*(?P<vx>{num}).*?"
            rf"VY\s*=\s*(?P<vy>{num}).*?"
            rf"VZ\s*=\s*(?P<vz>{num})",
            rec,
            re.DOTALL,
        )
        if m:
            vectors_europa.append([
                float(m.group("x")),
                float(m.group("y")),
                float(m.group("z")),
                float(m.group("vx")),
                float(m.group("vy")),
                float(m.group("vz")),
            ])

    print(vectors_europa[0])
    print(len(vectors_europa))

    # Ganymede query

    params_ganymede = {
        "format": "json",
        "COMMAND": "'503'", # Ganymede ID
        "MAKE_EPHEM": "'YES'",
        "OBJ_DATA": "'NO'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@599'", # Jupiter Center
        "START_TIME": "'1600-01-11 00:00:00'",
        "STOP_TIME": "'2200-01-09 00:00:00'",
        "STEP_SIZE": "'10 d'",
        "VEC_TABLE": "'2'",
        "VEC_CORR": "'NONE'",
        "OUT_UNITS": "'KM-S'",
        "TIME_TYPE": "'TDB'",
        "CAL_FORMAT": "'JD'",
        "CSV_FORMAT": "'NO'",
        "VEC_LABELS": "'YES'",
    }

    r_ganymede = requests.get(url, params=params_ganymede, timeout=30)
    r_ganymede.raise_for_status()
    data_ganymede = r_ganymede.json()

    if "error" in data_ganymede:
        raise RuntimeError(data_ganymede["error"])

    text_ganymede = data_ganymede["result"]

    block_ganymede = text_ganymede.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()

    records_ganymede = re.split(r"(?=\s*\d+\.\d+\s*=)", block_ganymede)

    vectors_ganymede = []
    for rec in records_ganymede:
        m = re.search(
            rf"(?P<jd>{num}).*?"
            rf"X\s*=\s*(?P<x>{num}).*?"
            rf"Y\s*=\s*(?P<y>{num}).*?"
            rf"Z\s*=\s*(?P<z>{num}).*?"
            rf"VX\s*=\s*(?P<vx>{num}).*?"
            rf"VY\s*=\s*(?P<vy>{num}).*?"
            rf"VZ\s*=\s*(?P<vz>{num})",
            rec,
            re.DOTALL,
        )
        if m:
            vectors_ganymede.append([
                float(m.group("x")),
                float(m.group("y")),
                float(m.group("z")),
                float(m.group("vx")),
                float(m.group("vy")),
                float(m.group("vz")),
            ])

    print(vectors_ganymede[0])
    print(len(vectors_ganymede))


    # Callisto query

    params_callisto = {
        "format": "json",
        "COMMAND": "'504'", # Callisto ID
        "MAKE_EPHEM": "'YES'",
        "OBJ_DATA": "'NO'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@599'", # Jupiter Center
        "START_TIME": "'1600-01-11 00:00:00'",
        "STOP_TIME": "'2200-01-09 00:00:00'",
        "STEP_SIZE": "'10 d'",
        "VEC_TABLE": "'2'",
        "VEC_CORR": "'NONE'",
        "OUT_UNITS": "'KM-S'",
        "TIME_TYPE": "'TDB'",
        "CAL_FORMAT": "'JD'",
        "CSV_FORMAT": "'NO'",
        "VEC_LABELS": "'YES'",
    }

    r_callisto = requests.get(url, params=params_callisto, timeout=30)
    r_callisto.raise_for_status()
    data_callisto = r_callisto.json()

    if "error" in data_callisto:
        raise RuntimeError(data_callisto["error"])

    text_callisto = data_callisto["result"]

    block_callisto = text_callisto.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()

    records_callisto = re.split(r"(?=\s*\d+\.\d+\s*=)", block_callisto)

    vectors_callisto = []
    for rec in records_callisto:
        m = re.search(
            rf"(?P<jd>{num}).*?"
            rf"X\s*=\s*(?P<x>{num}).*?"
            rf"Y\s*=\s*(?P<y>{num}).*?"
            rf"Z\s*=\s*(?P<z>{num}).*?"
            rf"VX\s*=\s*(?P<vx>{num}).*?"
            rf"VY\s*=\s*(?P<vy>{num}).*?"
            rf"VZ\s*=\s*(?P<vz>{num})",
            rec,
            re.DOTALL,
        )
        if m:
            vectors_callisto.append([
                float(m.group("x")),
                float(m.group("y")),
                float(m.group("z")),
                float(m.group("vx")),
                float(m.group("vy")),
                float(m.group("vz")),
            ])

    print(vectors_callisto[0])
    print(len(vectors_callisto))

    ephemerides = (vectors_io, vectors_europa, vectors_ganymede, vectors_callisto)

    return ephemerides