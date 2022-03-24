
materials = []

materials.append({'name': 'UDGRP', 'units': 'MPa-mm-Mg', 'rho': 2e-09,
                  'E1': 38000, 'E2': 8500, 'v12': 0.28, 'G12': 3400,
                  'XT': 1150, 'YT': 40, 'XC': 700, 'YC': 120, 'S12': 60, 'S23': 40, 'f12': -0.5})


def get(matname):
    for m in materials:
        if m['name'] == matname:
            return m
    return False
