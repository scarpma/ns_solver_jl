# defining units
U = 0.81
rho = 1060
d = 0.027

# derived units
rhoU2 = rho * U**2.
rhoU2d = rho * U**2. * d
t = d / U

# physical constants
nv = 4164
nt = 8229
tot_area_adim = 20.7704
tot_area = tot_area_adim * d**2.
area_adim = tot_area_adim / nt
area = area_adim * d**2.
E = 500e3
E_adim = E / rhoU2
h = 0.002
h_adim = h / d
E_2d = E * h
E_2d_adim = E_2d / rhoU2d
rho_body_adim = 1.5
rho_body = rho_body_adim * rho
m_tot = rho_body * h * tot_area
m_tot_adim = m_tot / (rho * d**3.)
m = m_tot / nv
m_adim = m / (rho * d**3.)
pois = 0.5
B = E * h**3. / (12 * (1 - pois**2.))
B_adim = B / (rhoU2 * d**3.)
kb = 2. * B / (3.**0.5)
kb_adim = kb / (rhoU2 * d**3.)
nu_crit = 2 * (m * E * h)**0.5
nu_crit_adim = nu_crit / (rho * U * d**2.)
reaction_time = 2 * m / nu_crit
reaction_time_adim = reaction_time / (d / U)




sepLen = 50
sep1 = "="*sepLen + "\n"
sep2 = "=" + "-"*(sepLen-2) + "="
i = 1
lineT = "{:1d}) {:6s} = {:<9.3g}   {:10s}"
lineT2 =    "  {:10s} = {:<9.3g}   {:10s}"
lineT3 = "  {:10s} = {:<9.3g}   {:10s}   {:<9.3g}"



print("dimension defining parameters (flow)")
print(sep1)

print(lineT.format(i, "U", U, "m s^-1"))
print(sep2)
i += 1

print(lineT.format(i, "rho", rho, "Kg m^-3"))
print(sep2)
i += 1

print(lineT.format(i, "d", d, "m"))
print(sep2)
i += 1

print()
print()
print("derived units")
print(sep1)

print(lineT2.format("t = d / U", t, "s"))
print(sep2)

print(lineT2.format("rho U^2", rhoU2, "Kg s^2 m^-1 = Pa"))
print(sep2)

print(lineT2.format("rho U^2 d", rhoU2d, "Kg s^2"))
print(sep2)


print()
print()
print("physical parameters")
print(sep1)

print(lineT3.format("nv", nv, "", nv))
print(sep2)

print(lineT3.format("nt", nt, "", nt))
print(sep2)

print(lineT3.format("tot_area", tot_area, "m^2", tot_area_adim))
print(sep2)

print(lineT3.format("area", area, "m^2", area_adim))
print(sep2)

print(lineT3.format("E", E, "Pa", E_adim))
print(sep2)

print(lineT3.format("h", h, "m", h_adim))
print(sep2)

print(lineT3.format("E_2d", E_2d, "Pa", E_2d_adim))
print(sep2)

print(lineT3.format("rho_body", rho_body, "kg m^-3", rho_body_adim))
print(sep2)

print(lineT3.format("pois ratio", pois, "", pois))
print(sep2)

print(lineT3.format("B", B, "Pa m^3", B_adim))
print(sep2)

print(lineT3.format("kb", kb, "Pa m^3", kb_adim))
print(sep2)

print(lineT3.format("m_tot", m_tot, "Kg", m_tot_adim))
print(sep2)

print(lineT3.format("m", m, "Kg", m_adim))
print(sep2)

print(lineT3.format("nu crit", nu_crit, "N s m^-1", nu_crit_adim))
print(sep2)

print(lineT3.format("react time", reaction_time, "s", reaction_time_adim))
print(sep2)


