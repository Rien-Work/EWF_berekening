import xlwings
from openpyxl import load_workbook
from CoolProp.HumidAirProp import HAPropsSI as lucht
import pvlib
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import time
import math

def weergegevens_ophalen():
    referentie_jaar = "2020"


    dateparse = lambda x: datetime.strptime(x, '%d-%m-%Y %H:%M:%S')
    dateparse1 = lambda x: datetime.strptime(x, '%d-%m-%Y %H:%M')

    # date_cols = ["DATE"]
    weer = pd.read_csv("weerdatapunten.csv", index_col='DATE', parse_dates=['DATE'], date_parser=dateparse)
    wind = pd.read_csv("weerdatapunten_v3.csv", index_col='DATE', parse_dates=['DATE'], date_parser=dateparse1)

    klimaatcascade = pd.read_csv("oud/klimaatcascade.csv")
    ventecdak = pd.read_csv("oud/ventecdak.csv")


    weer = pd.concat([weer, wind["wind_snelheid"], wind["zon_j_cm2"]], axis=1)

    periode = weer.loc[f'{referentie_jaar}-08-31 0:00:00':f'{referentie_jaar}-08-31 23:00:00']

    periode['dagvanweek'] = periode.index.dayofweek
    periode['uur'] = periode.index.hour

    # gegevens voor vakantie
    periode['jaar'] = periode.index.year
    periode['maand'] = periode.index.month
    periode['dag'] = periode.index.day
    lijst = periode.values.tolist()
    return lijst, periode

def bezetting_bepalen(periode):
    bezetting = []
    weekend = 0                 #bezetting in procenten in het weekend
    avond = 0.1                 #bezetting in de avond in procenten
    overdag = 0.5               #bezetting overdag
    df_bezetting = pd.DataFrame(periode, columns=['dagvanweek', 'uur'])
    #print(df_bezetting)
    lijst = df_bezetting.values.tolist()
    for x in lijst:
        if x[0] > 5:
            bezetting.append(weekend)
        elif 8 < x[1] < 18:
            bezetting.append(overdag)
        elif 17 < x[1] < 22:
            bezetting.append(avond)
        elif 21 < x[1] < 23:
            bezetting.append(0)
        else:
            bezetting.append(0)

    periode['bezetting'] = bezetting
    return periode

def debiet(periode):
    debiet_lijst = []
    debiet_m3_h = 122_000  # debiet in m3/h
    df_debiet = pd.DataFrame(periode, columns=['bezetting'])
    print(df_debiet)
    lijst = df_debiet.values.tolist()
    for x in lijst:
        debiet_bezetting = debiet_m3_h * x[0]
        debiet_lijst.append(debiet_bezetting)
    periode['lucht_debiet_m3_h'] = debiet_lijst
    return periode

def druk(periode):
    druk_lijst = []
    druk = 200  # debiet in m3/h
    df_druk = pd.DataFrame(periode, columns=['bezetting'])
    #print(df_debiet)
    lijst = df_druk.values.tolist()
    for x in lijst:
        druk_bezetting = druk * x[0]
        druk_lijst.append(druk_bezetting)
    periode['benodigde_druk_Pa'] = druk_lijst
    return periode

def lucht_dichtheid_kg_m3(periode):
    luchtdichtheid_lijst = []
    df_luchtdichtheid = pd.DataFrame(periode, columns=['temperatuur', 'RH'])
    lijst = df_luchtdichtheid.values.tolist()
    for x in lijst:
        kg_m3 = (1 / lucht('Vda', 'T', x[0] + 273.15, 'P', 101325, 'R', x[1]/100))
        luchtdichtheid_lijst.append(round(kg_m3,3))
    periode['luchtdichtheid_kg_m3'] = luchtdichtheid_lijst
    return periode

def kg_kg_vocht_in_lucht_berekenen(periode):
    luchtvochtigheid_lijst = []
    df_luchtvocht = pd.DataFrame(periode, columns=['temperatuur', 'RH'])
    lijst = df_luchtvocht.values.tolist()
    for x in lijst:
        w_lucht = lucht('W', 'T', x[0] + 273.15, 'P', 101325, 'R', x[1]/100)
        luchtvochtigheid_lijst.append(round(w_lucht,4))
    periode['luchtvocht_kg_kg'] = luchtvochtigheid_lijst
    return periode

def berekenen_klimaatcascade(invoer_t_air_in, invoer_m_air_in_m3_h, invoer_rho_air, invoer_x_air_in_g_kg, invoer_r_water_air,invoer_t_wat_in, invoer_t_air_out_side):
    #voorbewerken waardes
    invoer_x_air_in_g_kg = invoer_x_air_in_g_kg*1000

    e = 2.71828182845904
    klimaatcascade = pd.DataFrame()

    # specificaties sproeiers
    d_wat = 0.581  # mm
    d_vmd = 1.708  # mm
    d_smd = 1.377  # mm
    w_vt = 9.53  # l/min

    # klimaatcascade
    height = 35  # meter
    height_total = height

    width = 1.30  # meter
    length = 1.3  # meter

    dh = height / 200
    time = 0
    position = 0

    # air
    rho_air = invoer_rho_air  # kg/m3
    m_air_in_m3_h = invoer_m_air_in_m3_h  # m3/h
    m_air_in_kg_s = m_air_in_m3_h * rho_air / 3600  # kg/s

    # heat tranfer
    fall_velocity = 1.741 * d_wat + 0.1623  # m/s
    h_c = (2 + 0.6 * 0.71 ** (2 / 3) * math.sqrt(
        fall_velocity * d_wat / 1000 / 0.000015)) * 0.025 / d_wat * 1000  # w/m2K
    h_w = 3.98 * ((m_air_in_kg_s / 1.2) / width / length) ** 0.8  # w/m2K


    t_air_in = invoer_t_air_in  # graC
    t_air_in_begin = t_air_in
    t_air_out_side = invoer_t_air_out_side

    x_sat_air_in = 7.22  # g/kg

    c_air = 1012  # j/kgk
    v_air = (m_air_in_m3_h / 3600) / (width * length)  # m/s

    # vapour
    x_air_in_g_kg = invoer_x_air_in_g_kg
    x_air_in = x_air_in_g_kg / 1000  # kg/kg

    rho_v_in = rho_air * x_air_in  # kg/m3
    c_vap = 2020  # j/kg/k

    L_wat = 2257000  # J/kg

    # water
    r_water_air = invoer_r_water_air  # kg/kg
    m_wat_in = r_water_air * m_air_in_kg_s  # kg water per kg lucht

    t_wat_in = invoer_t_wat_in  # graC

    rho_a_wat_in_kg_m3 = 0.1804714474  # kg/m3
    c_wat = 4184  # j/kgK
    v_wat = fall_velocity + v_air

    # eerste variable berekenen voor eerste rij
    p_vapour_in_air_Pa = rho_v_in * 462 * (t_air_in + 273.15)
    p_sat_water_Pa = 100 * e ** (18.956 - (4030.18) / (t_wat_in + 235))

    p_vapour_in_air_Pa_power = rho_v_in * 462 * (t_air_in + 273.15)

    # heat transfer cascade
    number_of_droplets_per_m3 = rho_a_wat_in_kg_m3 / ((1000 * 4 / 3) * math.pi * (((d_vmd / 2) / 1000) ** 3))
    surface_area_droplets_per_m3 = number_of_droplets_per_m3 * 4 * math.pi * (d_smd / 2 / 1000) ** 2
    mass_tranf_per_m3_kg_m3s = 0.0000062 * 0.00089 * h_c * surface_area_droplets_per_m3 * (
                p_vapour_in_air_Pa - p_sat_water_Pa)
    P_latent_droplets_W_m3 = mass_tranf_per_m3_kg_m3s * (L_wat + c_wat * (100 - t_wat_in) - c_vap * (100 - t_air_in))
    P_sensible_droplets_W_m3 = (t_air_in - t_wat_in) * surface_area_droplets_per_m3 * h_c

    # heat transfer walls
    surface_area_per_m3 = (2 * width + 2 * length) / (width * length)
    mass_tranf_wall_per_m3_kg_m3 = 0.0000062 * 0.00089 * h_w * surface_area_per_m3 * (
                p_vapour_in_air_Pa - p_sat_water_Pa)
    P_latent_walls_W_m3 = mass_tranf_wall_per_m3_kg_m3 * (L_wat + c_wat * (100 - t_wat_in) - c_vap * (100 - t_air_in))
    P_sensible_walls_W_m3 = (t_air_in - t_wat_in) * surface_area_per_m3 * h_w

    # total heat transfer
    P_latten_total_W_m3 = P_latent_droplets_W_m3 + P_latent_walls_W_m3
    P_sensible_total_W_m3 = P_sensible_walls_W_m3 + P_sensible_droplets_W_m3

    # lijsten
    # air_vapour
    height_in_cascade = [height]
    time_list = [time]
    temp_air_c_list = [t_air_in]
    rho_vapour_kg_m3_list = [rho_v_in]
    p_vapour_in_air_Pa_list = [p_vapour_in_air_Pa]

    # water
    rho_a_wat_in_kg_m3_list = [rho_a_wat_in_kg_m3]
    temp_water_C_list = [t_wat_in]
    p_sat_water_Pa_list = [p_sat_water_Pa]

    # heat tranfer cascade
    number_of_droplets_per_m3_list = [number_of_droplets_per_m3]
    surface_area_droplets_per_m3_list = [surface_area_droplets_per_m3]
    mass_tranf_per_m3_kg_m3s_list = [mass_tranf_per_m3_kg_m3s]
    P_latent_droplets_W_m3_list = [P_latent_droplets_W_m3]
    P_sensible_droplets_W_m3_list = [P_sensible_droplets_W_m3]

    # heat transfer walls
    surface_area_per_m3_list = [surface_area_per_m3]
    mass_tranf_wall_per_m3_kg_m3_list = [mass_tranf_wall_per_m3_kg_m3]
    P_latent_walls_W_m3_list = [P_latent_walls_W_m3]
    P_sensible_walls_W_m3_list = [P_sensible_walls_W_m3]

    # total heat transfer
    P_latten_total_W_m3_list = [P_latten_total_W_m3]
    P_sensible_total_W_m3_list = [P_sensible_total_W_m3]

    while height > 0:
        height = height - dh
        time = time + dh / v_air
        # air+vapour
        rho_v_in = rho_vapour_kg_m3_list[position] - dh / v_air * mass_tranf_per_m3_kg_m3s_list[position]

        t_air_in = 1 / (c_vap * rho_v_in + c_air * rho_air) * (
                    c_vap * rho_vapour_kg_m3_list[position] * temp_air_c_list[position] + c_air * rho_air *
                    temp_air_c_list[position] + L_wat * (rho_vapour_kg_m3_list[position] - rho_v_in) - dh / v_air * (
                                P_latten_total_W_m3_list[position] + P_sensible_total_W_m3_list[position]))
        p_vapour_in_air_Pa = rho_v_in * 462 * (t_air_in + 273.15)

        # water
        rho_a_wat_in_kg_m3 = rho_a_wat_in_kg_m3_list[position] + dh / v_wat * mass_tranf_per_m3_kg_m3s_list[position]
        t_wat_in = 1 / (c_wat * rho_a_wat_in_kg_m3) * (
                    c_wat * rho_a_wat_in_kg_m3_list[position] * temp_water_C_list[position] + dh / v_wat * (
                        P_latten_total_W_m3_list[position] + P_sensible_total_W_m3_list[position]))
        p_sat_water_Pa = 100 * e ** (18.956 - (4030.18) / (t_wat_in + 235))

        # heat tranfer cacade
        number_of_droplets_per_m3 = rho_a_wat_in_kg_m3 / ((1000 * 4 / 3) * math.pi * (((d_vmd / 2) / 1000) ** 3))
        surface_area_droplets_per_m3 = number_of_droplets_per_m3 * 4 * math.pi * (d_smd / 2 / 1000) ** 2
        mass_tranf_per_m3_kg_m3s = 0.0000062 * 0.00089 * h_c * surface_area_droplets_per_m3 * (
                    p_vapour_in_air_Pa - p_sat_water_Pa)
        P_latent_droplets_W_m3 = mass_tranf_per_m3_kg_m3s * (
                    L_wat + c_wat * (100 - t_wat_in) - c_vap * (100 - t_air_in))
        P_sensible_droplets_W_m3 = (t_air_in - t_wat_in) * surface_area_droplets_per_m3 * h_c
        # heat transfer walls
        surface_area_per_m3 = (2 * width + 2 * length) / (width * length)
        mass_tranf_wall_per_m3_kg_m3 = 0.0000062 * 0.00089 * h_w * surface_area_per_m3 * (
                    p_vapour_in_air_Pa - p_sat_water_Pa)
        P_latent_walls_W_m3 = mass_tranf_wall_per_m3_kg_m3 * (
                    L_wat + c_wat * (100 - t_wat_in) - c_vap * (100 - t_air_in))
        P_sensible_walls_W_m3 = (t_air_in - t_wat_in) * surface_area_per_m3 * h_w

        # total heat transfer
        P_latten_total_W_m3 = P_latent_droplets_W_m3 + P_latent_walls_W_m3
        P_sensible_total_W_m3 = P_sensible_walls_W_m3 + P_sensible_droplets_W_m3

        # toevoegen aan lijst
        height_in_cascade.append(height)
        time_list.append(round(time, 3))
        rho_vapour_kg_m3_list.append(rho_v_in)
        temp_air_c_list.append(t_air_in)
        p_vapour_in_air_Pa_list.append(p_vapour_in_air_Pa)
        # water
        rho_a_wat_in_kg_m3_list.append(rho_a_wat_in_kg_m3)
        temp_water_C_list.append(round(t_wat_in, 3))
        p_sat_water_Pa_list.append(p_sat_water_Pa)

        # heat tranfer cacade
        number_of_droplets_per_m3_list.append(number_of_droplets_per_m3)
        surface_area_droplets_per_m3_list.append(surface_area_droplets_per_m3)
        mass_tranf_per_m3_kg_m3s_list.append(mass_tranf_per_m3_kg_m3s)
        P_latent_droplets_W_m3_list.append(P_latent_droplets_W_m3)
        P_sensible_droplets_W_m3_list.append(P_sensible_droplets_W_m3)

        # heat transfer walls
        surface_area_per_m3_list.append(surface_area_per_m3)
        mass_tranf_wall_per_m3_kg_m3_list.append(mass_tranf_wall_per_m3_kg_m3)
        P_latent_walls_W_m3_list.append(P_latent_walls_W_m3)
        P_sensible_walls_W_m3_list.append(P_sensible_walls_W_m3)

        P_latten_total_W_m3_list.append(P_latten_total_W_m3)
        P_sensible_total_W_m3_list.append((P_sensible_total_W_m3))

        # print(position)

        position = position + 1  # om positie te bepalen en juiste slectie te maken in lijst

    # gegevens toevoegen aan dataframe
    # air+vapour
    klimaatcascade['time_in_s'] = time_list
    klimaatcascade['height_in_cascade_m'] = height_in_cascade
    klimaatcascade['temp_air_C'] = temp_air_c_list
    klimaatcascade['rho_vouper_kg_m3'] = rho_vapour_kg_m3_list

    # water
    klimaatcascade['rho_a_wat)kg_m3'] = rho_a_wat_in_kg_m3_list
    klimaatcascade['temp_water_c'] = temp_water_C_list
    klimaatcascade['P_sat_water_Pa'] = p_sat_water_Pa_list


    temp_air_out = float(klimaatcascade['temp_air_C'].tail(1))  # temperatuur lucht aan voet van climaatcascade
    temp_water_out = float(klimaatcascade['temp_water_c'].tail(1))
    vouper_air_out = float(klimaatcascade['P_sat_water_Pa'].tail(1))

    pd_hydraulisch_druk = height_total * m_wat_in * 9.81 / ((width * length) * (w_vt + v_air))
    t_Kc_gem = (temp_air_out * 2 / 3) + ((t_air_in_begin + temp_air_out) / 2) * 1 / 3
    pd_thermische_druk = 1.293 * 9.81 * height_total * (((273 / (273 + t_air_out_side)) - ((273 / (273 + t_Kc_gem)))))
    pd_druk_totaal = pd_hydraulisch_druk - pd_thermische_druk

    RH_air_out = vouper_air_out / (100 * e**(18.956 - (4030.18) / (temp_air_out + 235)))

    total_cooling_power_sensible =(t_air_in_begin - temp_water_out) * m_air_in_kg_s * c_air / 1000
    total_cooling_power_latent = (p_vapour_in_air_Pa_power-vouper_air_out)*m_air_in_kg_s/rho_air*L_wat/1000
    total_cooling_power = total_cooling_power_sensible + total_cooling_power_latent

    if RH_air_out > 1:
        RH_air_out = 0.99

    klimaatcascade.to_csv("C:/Users/MSPE/PycharmProjects/EWF concept/klimaatcascade.csv")
    return temp_air_out, RH_air_out, pd_druk_totaal, temp_water_out, total_cooling_power, m_wat_in


def klimaatcascade_ventilator_verbruik(dp, debiet_m3_h, druk_instalatie):
    ventilator_rendement = 0.6                          #ventilator van de ventilatoren boven klimaatcascade
    benodigde_druk = druk_instalatie - dp               #instalatiedruk - druk opgewekt in klimaatcascade
    opgenomen_kw_ven = (debiet_m3_h*benodigde_druk)/(3600*1000*(ventilator_rendement))

    return opgenomen_kw_ven

def klimaatcascade_lucht_voorverwarmen_twincoil(t_lucht_buiten, RH_lucht_buiten, debiet_m3_h):

    t_lucht_buiten = ((21 - t_lucht_buiten)*0.2) + t_lucht_buiten

    lucht_voorverwarmen_tot = 3                                                                                     #lucht wordt voorverwarmd tot 3 graden om bevriezen te voorkomen
    RH_uit = 0
    if t_lucht_buiten < lucht_voorverwarmen_tot:
        t_lucht_buiten = t_lucht_buiten + 273.15

        w_in = lucht('W', 'T', t_lucht_buiten, 'P', 101325, 'R', RH_lucht_buiten/100)
        h_in = lucht('H', 'T', t_lucht_buiten, 'P', 101300, 'W', w_in)
        h_uit = lucht('H', 'T', lucht_voorverwarmen_tot + 273.15, 'P', 101300, 'W', w_in)
        RH_uit = lucht('R', 'T', t_lucht_buiten, 'P', 101325, 'W', w_in)

        volumestroom = debiet_m3_h / 3600                                                                           #van m3/h naar m3/s
        h_voorverwarmen_entalpie = h_in - h_uit                                                                     #verschil in enthalpie van luchtstroom
        kg_m3 = (1 / lucht('Vda', 'T', t_lucht_buiten, 'P', 101325, 'R', RH_lucht_buiten/100))                          #massa van lucht berekenen m3/kg/1=kg/m3
        kW_voorverwarmen = ((kg_m3 * h_voorverwarmen_entalpie * volumestroom) / 1000) * -1                          #kg/m3 * kJ/kg * m3/s = kJ/s
        lucht_voorverwarmen_tot = 3
    else:                                                                                                           #temp onder voorverwarm voorwaarde, nul verbruik
        kW_voorverwarmen = 0
        lucht_voorverwarmen_tot = t_lucht_buiten
    return kW_voorverwarmen, RH_uit, lucht_voorverwarmen_tot

def klimaatcascade_lucht_naverwarmen(t_lucht_uit, RH_lucht_uit, debiet_m3_h):
    lucht_naverwarmen_tot = 18                                                                                      #lucht wordt voorverwarmd tot 3 graden om bevriezen te voorkomen
    lucht_nakoelen_tot = 22

    if RH_lucht_uit > 99:                                                                                           #uit excel van Ben komt wel eens een RH van meer dan 100%, dat wodt hier terug geschoefd
        RH_lucht_uit = 99

    if t_lucht_uit < lucht_naverwarmen_tot:                                                                         #bepalen of lucht opgewarmd moet worden of niet
        t_lucht_uit = t_lucht_uit + 273.15

        w_in = lucht('W', 'T', t_lucht_uit, 'P', 101325, 'R', RH_lucht_uit/100)                                         #kg water per kg lucht bepalen
        h_in = lucht('H', 'T', t_lucht_uit, 'P', 101300, 'W', w_in)                                                 #entalpie van inkomende lucht bepalen
        h_uit = lucht('H', 'T', lucht_naverwarmen_tot + 273.15, 'P', 101300, 'W', w_in)                             #entalpie van uitlkomende lucht bepalen
        RH_uit = lucht('R', 'T', lucht_naverwarmen_tot + 273.15, 'P', 101325, 'W', w_in)                                     #luchtvochtigheid uitgaande lucht
        RH_uit_na_voorwerwarmeer = RH_uit
        volumestroom = debiet_m3_h / 3600                                                                           #van m3/h naar m3/s
        h_voorverwarmen_entalpie = h_in - h_uit                                                                     #verschil in enthalpie van luchtstroom
        kg_m3 = (1 / lucht('Vda', 'T', t_lucht_uit, 'P', 101325, 'R', RH_lucht_uit/100))                                #massa van lucht berekenen m3/kg/1=kg/m3
        kW_naverwarmen = ((kg_m3 * h_voorverwarmen_entalpie * volumestroom) / 1000) * -1                          #kg/m3 * kJ/kg * m3/s = kJ/s
        kW_nakoelen = 0

    elif t_lucht_uit > lucht_nakoelen_tot:                                                                         #bepalen of lucht opgewarmd moet worden of niet
        t_lucht_uit = t_lucht_uit + 273.15

        w_in = lucht('W', 'T', t_lucht_uit, 'P', 101325, 'R', RH_lucht_uit/100)                                         #kg water per kg lucht bepalen
        h_in = lucht('H', 'T', t_lucht_uit, 'P', 101300, 'W', w_in)                                                 #entalpie van inkomende lucht bepalen
        h_uit = lucht('H', 'T', lucht_nakoelen_tot + 273.15, 'P', 101300, 'W', w_in)                             #entalpie van uitlkomende lucht bepalen
        RH_uit = lucht('R', 'T', lucht_nakoelen_tot, 'P', 101325, 'W', w_in)                                     #luchtvochtigheid uitgaande lucht
        RH_uit_na_voorwerwarmeer = RH_uit * 10000
        volumestroom = debiet_m3_h / 3600                                                                           #van m3/h naar m3/s
        h_voorverwarmen_entalpie = h_in - h_uit                                                                     #verschil in enthalpie van luchtstroom
        kg_m3 = (1 / lucht('Vda', 'T', t_lucht_uit, 'P', 101325, 'R', RH_lucht_uit/100))                                #massa van lucht berekenen m3/kg/1=kg/m3
        kW_nakoelen = ((kg_m3 * h_voorverwarmen_entalpie * volumestroom) / 1000)
        kw_naverwarmen = 0

    else:                                                                                                           #temp onder voorverwarm voorwaarde, nul verbruik
        kW_naverwarmen = 0                                                                                        #als lucht 18 graC wordt er niet naverwarmt
        kW_nakoelen = 0

    return kW_naverwarmen, RH_uit_na_voorwerwarmeer, kW_nakoelen

def klimaatcascade_pomp(m_wat_in):
    return (m_wat_in*3600) / 60 * (3 + 0.6) / 540       #pomp energie in kW

def klimaatcascade_water_verwarmen(t_wat_out, m_wat_in, t_wat_set):
    if t_wat_out < t_wat_set:       #nodig om water te verwarmen tot setpoint water temp in waterbak
        t_wat_verschil = t_wat_set - t_wat_out
        return t_wat_verschil * 4.2 * m_wat_in #delta T [C] * 4,3 kj/kgK * massastroom [kg/s]

def opbrengs_pvt(solar_radiation):
    surface_pvt = 1.8 * 100                                         #honderd panelen van 1.8 m2
    rendement = 0.7
    totale_radiatie = ((solar_radiation* 100**2 * surface_pvt)/1000)/3600  #instraling in kJ/s
    opbrenngst_pvt_totaal = totale_radiatie * rendement
    return opbrenngst_pvt_totaal                            #kj/s

def twin_coil(debiet_m3_h, t_lucht_binnen, t_lucht_buiten, RH_buiten):
    rendement_twin_coil = 0.2                                   #rendement is 20 procent
    t_binnen = t_lucht_binnen + 273.15                                #binnen temperatuur uitgaande van 21 procent
    w_in = lucht('W', 'T', t_lucht_binnen, 'P', 101325, 'R', 0.6)  # kg water per kg lucht bepalen
    h_in = lucht('H', 'T', t_lucht_binnen, 'P', 101300, 'W', w_in)  # entalpie van inkomende lucht bepalen
    h_uit = lucht('H', 'T', t_lucht_buiten + 273.15, 'P', 101300, 'R',w_in)  # entalpie van uitlkomende lucht bepalen
    volumestroom = debiet_m3_h / 3600  # van m3/h naar m3/s
    h_voorverwarmen_entalpie = h_in - h_uit  # verschil in enthalpie van luchtstroom
    kg_m3 = (1 / lucht('Vda', 'T', t_lucht_binnen, 'P', 101325, 'R', 0.6))  # massa van lucht berekenen m3/kg/1=kg/m3
    kW_twin_coil = ((kg_m3 * h_voorverwarmen_entalpie * volumestroom) / 1000) * -1
    return kW_twin_coil * rendement_twin_coil

def beloke_schoorsteen(opbrengt_pvt, t_lucht_buiten, t_lucht_binnen, debiet_m3_h):
    w_binnen = lucht('W', 'T', t_lucht_binnen+273.15, 'P', 101325, 'R', 0.6)  # kg water per kg lucht bepalen
    h_binnen = lucht('H', 'T', t_lucht_binnen+273.15, 'P', 101300, 'W', w_binnen)  # entalpie van inkomende lucht bepalen
    kg_m3 = (1 / lucht('Vda', 'T', t_lucht_binnen+273.15, 'P', 101325, 'R', 0.6))
    volumestroom = debiet_m3_h / 3600 #m3/s
    kg_s = volumestroom * kg_m3
    vermoge_lucht = kg_s * h_binnen + (opbrengt_pvt *1000)
    h_lucht_na = vermoge_lucht/kg_s
    t_lucht_schoorsteen = lucht('T', 'H', h_lucht_na, 'P', 101300, 'W',w_binnen) - 273.15 # entalpie van uitlkomende lucht bepalen
    Pa_druk_schoorsteen = 3465 * 20 * ((1/(t_lucht_buiten+273.15)) - (1/(t_lucht_schoorsteen+273.15)))
    return Pa_druk_schoorsteen

def ventec_dak(wind):
    ventecdak = pd.read_csv("oud/ventecdak.csv")
    if wind > 8:
        wind = 8
    positie = int(ventecdak[ventecdak["windsnelheid"] == wind].index.values)
    wind_druk = ventecdak.iloc[positie, 1]
    return wind_druk

def berekening_klimaatcascade(periode):
    #gevens uit df halen voor waterhoeveelheden klimaatcascade
    klimaatcascade = pd.read_csv("klimaatcascade_water_hoeveelheid.csv")

    df_klimaatcascade = pd.DataFrame(periode, columns=['temperatuur', 'RH', 'luchtdichtheid_kg_m3', "lucht_debiet_m3_h", "luchtvocht_kg_kg"])
    klimaatcascade_lijst = []
    t_air_uit_list = []
    RH_air_out_list = []
    dp_list = []
    t_wat_out_list = []
    total_cooling_power_list = []
    m_wat_in_list = []
    lucht_voorverwarmen_tot_list = []
    #voorverwarmen van lucht:
    kW_voorverwarmen_lucht_list = []
    RH_uit_navoorverwarmer_list = []
    #naverwarmen
    kW_naverwarmer_list = []
    kW_nakoeler_list = []
    RH_inblaas_gebouw_list = []

    #pomp energie
    kW_pomp_list = []
    voortgang = 0                                           #variable voor percentage uit te rekenen van berekening
    lijst = df_klimaatcascade.values.tolist()
    count = 0
    for x in lijst:
        count += 1
        print(str(round((count/len(lijst))*100, 3)) + "%") #bereken op hoeveel procent de berekening is
        #print(x[3])
        if x[3] == 0:                                      #als luchtdebiet 0 wordt deze periode niet berekend
            #print("luchtdebiet:")
            #print(x[3])
            t_air_uit_list.append(0)
            RH_air_out_list.append(0)
            dp_list.append(0)
            t_wat_out_list.append(0)
            total_cooling_power_list.append(0)
            m_wat_in_list.append(0)

            kW_naverwarmer_list.append(0)
            kW_nakoeler_list.append(0)
            RH_inblaas_gebouw_list.append(0)

            kW_voorverwarmen_lucht_list.append(0)
            RH_uit_navoorverwarmer_list.append(0)
            lucht_voorverwarmen_tot_list.append(0)

            kW_pomp_list.append(0)
        elif x[3] > 0:                                  #als luchtdebiet groter is dan 0 wordt deze doorgerekend
            #water hoeveelheid selecteren
            temp_select = round(x[0])
            positie = int(klimaatcascade[klimaatcascade["na_twin"] == temp_select].index.values)
            water_in = klimaatcascade.iloc[positie, 1]



            kW_voorverwarmen, RH_uit_na_naverwarmer_, lucht_voorverwarmen_tot = klimaatcascade_lucht_voorverwarmen_twincoil(x[0], x[1], x[3])
            t_air_uit, RH_air_out, dp, t_wat_out, total_cooling_power, m_wat_in = berekenen_klimaatcascade(lucht_voorverwarmen_tot, x[3], x[2], x[4], water_in, 13, 10)
            kW_naverwarmen, RH_uit_na_voorwerwarmeer_, kW_nakoelen = klimaatcascade_lucht_naverwarmen(t_air_uit, RH_air_out, x[3])

            kW_pomp_list.append(klimaatcascade_pomp(m_wat_in))  # pomp energie berekenen
            #print(m_wat_in)
            kW_naverwarmer_list.append(kW_naverwarmen)
            kW_nakoeler_list.append(kW_nakoelen)
            RH_inblaas_gebouw_list.append(RH_uit_na_voorwerwarmeer_*100)

            t_air_uit_list.append(t_air_uit)
            RH_air_out_list.append(RH_air_out)
            dp_list.append(dp)
            t_wat_out_list.append(t_wat_out)
            total_cooling_power_list.append(total_cooling_power)
            m_wat_in_list.append(m_wat_in)

            lucht_voorverwarmen_tot_list.append(lucht_voorverwarmen_tot)

            #gegevens van voorwarwarmen toevoegen aan lijst:
            kW_voorverwarmen_lucht_list.append(kW_voorverwarmen)
            RH_uit_navoorverwarmer_list.append(RH_uit_na_naverwarmer_)


            #print(t_air_uit)
    #gegevens toevoegen aan uitkomsten bestand
    periode['temperatuur_voor_klimaatcascade'] = lucht_voorverwarmen_tot_list
    periode['temperatuur_na_klimaatcascade'] = t_air_uit_list
    periode['RH_na_klimaatcascade'] = RH_air_out_list
    periode['totaal_drukverschil_na_klimaatcascade'] = dp_list
    periode['water_temp_na_klimaatcascade'] = t_wat_out_list
    periode['Totaal_koelvermogen_klimaatcascade'] = total_cooling_power_list

    #pomp energie bereken
    periode['pomp_energie_in_kw'] = kW_pomp_list


    #naverwarmen toevoegen aan uitkomsten bestand
    periode['naverwarmen_kW'] = kW_naverwarmer_list
    periode['nakoelen_kW'] = kW_nakoeler_list
    periode['RH_inblaas_gebouw'] = RH_inblaas_gebouw_list
    return periode


def berekenen_beloke_schoorsteen_ventec_dak(periode):
    druk_schoorsteen = []
    ventilator_kw = []
    df_beloke_schoorsteen = pd.DataFrame(periode, columns=['temperatuur', 'RH', 'zon_j_cm2', "lucht_debiet_m3_h", "wind_snelheid", 'benodigde_druk_Pa'])
    lijst = df_beloke_schoorsteen.values.tolist()
    for x in lijst:
        if x[3] > 0:
            opbrengst_pvt_ = opbrengs_pvt(x[2])
            druk_schoorsteen_ = beloke_schoorsteen(opbrengst_pvt_, x[0], 21, x[3])
            druk_schoorsteen.append(round(druk_schoorsteen_, 3))
            druk_ventec_dak = ventec_dak(x[4])
            druk_totaal = druk_schoorsteen_ + druk_ventec_dak
            if druk_totaal < x[5]:
                ventilator_kw_ = klimaatcascade_ventilator_verbruik(druk_totaal, x[3], x[5])
                ventilator_kw.append(round(ventilator_kw_, 3))
        if x[3] == 0:
            druk_schoorsteen.append(0)
            ventilator_kw.append(0)
    periode['luchtdruk_belokeschoorsteen_Pa'] = druk_schoorsteen
    periode['Ventilator_energie_schoorsteen_kW'] = ventilator_kw
    return periode
#output: luchtdruk_belokeschoorsteen_Pa, Ventilator_energie_schoorsteen_kW


#t_air_uit, RH_air_out, dp, t_wat_out, total_cooling_power, m_wat_in = klimaatcascade_excel_invoer(10,10_000,1.2,0.006,0.1680,13)

lijst,periode = weergegevens_ophalen()

periode = bezetting_bepalen(periode)

periode = debiet(periode)
periode = druk(periode)
periode = lucht_dichtheid_kg_m3(periode)
periode = kg_kg_vocht_in_lucht_berekenen(periode)
periode = berekening_klimaatcascade(periode)
periode = berekenen_beloke_schoorsteen_ventec_dak(periode)

periode.to_csv("C:/Users/MSPE/PycharmProjects/EWF concept/uitkomsten.csv")