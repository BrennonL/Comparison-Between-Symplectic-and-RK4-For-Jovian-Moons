import requests
import re

url = "https://ssd.jpl.nasa.gov/api/horizons.api"
num = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:E[+-]?\d+)?"

def extract_mu_km3_s2(text: str) -> float:
    m = re.search(
        r"GM\s*\(km\^3/s\^2\)\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?)",
        text
    )
    if not m:
        raise RuntimeError("Could not find GM in Horizons response.")
    return float(m.group(1))

def get_body_data(body_id: str):
    params = {
        "format": "json",
        "COMMAND": f"'{body_id}'",
        "MAKE_EPHEM": "'YES'",
        "OBJ_DATA": "'YES'",
        "EPHEM_TYPE": "'VECTORS'",
        "CENTER": "'500@599'",
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

    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()

    if "error" in data:
        raise RuntimeError(data["error"])

    text = data["result"]
    mu = extract_mu_km3_s2(text)

    block = text.split("$$SOE", 1)[1].split("$$EOE", 1)[0].strip()
    records = re.split(r"(?=\s*\d+\.\d+\s*=)", block)

    julian_dates = []
    vectors = []
    for rec in records:
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
            julian_dates.append(float(m.group("jd")))
            vectors.append([
                float(m.group("x")),
                float(m.group("y")),
                float(m.group("z")),
                float(m.group("vx")),
                float(m.group("vy")),
                float(m.group("vz")),
            ])

    return julian_dates, vectors, mu

def get_body_mu(body_id: str) -> float:
    params = {
        "format": "json",
        "COMMAND": f"'{body_id}'",
        "MAKE_EPHEM": "'NO'",
        "OBJ_DATA": "'YES'",
    }

    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()

    if "error" in data:
        raise RuntimeError(data["error"])

    return extract_mu_km3_s2(data["result"])

def get_reference_ephemerides():
    moon_ids = ("501", "502", "503", "504")

    moon_data = [get_body_data(body_id) for body_id in moon_ids]

    time_vec = moon_data[0][0]
    body_major_ephemerides = [vectors for _, vectors, _ in moon_data]
    ephemerides = [
        [body_vectors[step_index] for body_vectors in body_major_ephemerides]
        for step_index in range(len(time_vec))
    ]
    moon_mus = tuple(mu for _, _, mu in moon_data)

    Jupter_mu = get_body_mu("599")

    return ephemerides, time_vec, moon_mus, Jupter_mu
