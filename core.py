# Script by E. J. Ridley
import os
import numpy as np
import requests
import json
import warnings
from matplotlib import pyplot as plt, gridspec as gspec
from astropy.time import Time, TimeDelta
from astroplan import FixedTarget, Observer
from astroplan.plots import plot_altitude
from astroplan.moon import moon_illumination
from astropy.coordinates import SkyCoord, get_moon
import astropy.units as u


# ignore annoying plotting warnings
warnings.filterwarnings('ignore', 'linestyle')


plt.rcParams['figure.figsize'] = 12, 7
plt.rcParams['font.size'] = 10
plt.rcParams['axes.prop_cycle'] = plt.cycler(
    color=['#583ea3', '#ff7f00', '#4daf4a', '#f781bf', '#164608', '#377eb8', '#999999', '#ff1a1c']
)  # modified from https://gist.github.com/thriveth/8560036


# index relating marshall priorities to an internal rank system
RANK_INDEX = {
    'CRITICAL': 1,
    'HIGH': 2,
    'IMPORTANT': 3,
    'MEDIUM': 4,
    'USEFUL': 5,
    'LOW': 6
}


# reversed dictionary for convenience
REV_RANK_INDEX = {V: K for K, V in RANK_INDEX.items()}


# current date, if desired replace DATE with date for planning e.g. '2021-08-07'
DATE = str(Time.now().datetime.date())


# ePESSTO+ observation blocks
CLASSIFICATION_OBS = [
    {'mag_constraint': lambda m: m < 13, 'exp_time_s': 40, 'exec_time_s': 433, 'description': '< 13'},
    {'mag_constraint': lambda m: 13 <= m < 14.5, 'exp_time_s': 120, 'exec_time_s': 517, 'description': '13-14.5'},
    {'mag_constraint': lambda m: 14.5 <= m < 16, 'exp_time_s': 180, 'exec_time_s': 581, 'description': '14.5-16'},
    {'mag_constraint': lambda m: 16 <= m < 17.5, 'exp_time_s': 300, 'exec_time_s': 711, 'description': '16-17.5'},
    {'mag_constraint': lambda m: 17.5 <= m < 18.5, 'exp_time_s': 600, 'exec_time_s': 1031, 'description': '17.5-18.5'},
    {'mag_constraint': lambda m: 18.5 <= m < 19.5, 'exp_time_s': 900, 'exec_time_s': 1351, 'description': '18.5-19.5'},
    {'mag_constraint': lambda m: 19.5 <= m, 'exp_time_s': 1500, 'exec_time_s': 1971, 'description': '> 19.5'}
]


def plot_wrap(fig, *args, show=False, filename=None, transparent=False):
    if transparent:
        fig.patch.set_alpha(0)
        for ax in fig.axes:
            ax.patch.set_facecolor('white')
            ax.patch.set_alpha(1.)

    if filename is not None:
        fig.savefig(f'{filename}.pdf')

    if show:
        plt.show()

    plt.close(fig)

    return args


def make_directory(path):
    if not os.path.exists(path):
        try:
            os.mkdir(path)
        except OSError:
            print(f'Failed to create directory {path}.')


def strip_name(name):
    if name[:2] == 'AT':
        name = name[4:]
    elif name[:2] == 'SN':
        name = name[4:]
    elif name[:3] == 'ZTF':
        name = name[3:]

    return name


def observation_timings(site, date):
    midday = Time(f'{date} 12:00')
    sunset = site.sun_set_time(midday, which='next')
    sunrise = site.sun_rise_time(midday, which='next')
    evening_twilight = site.twilight_evening_astronomical(midday, which='next')
    morning_twilight = site.twilight_morning_astronomical(midday, which='next')

    return sunset, sunrise, evening_twilight, morning_twilight


def target_altitudes(targets, site, date, min_rank=1, max_rank=6):
    fig = plt.figure()
    gs = gspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])

    n_colors = sum([1 for d in plt.rcParams['axes.prop_cycle']])

    # get timings
    hour = TimeDelta(3600.0 * u.s)
    sunset, sunrise, evening_twilight, morning_twilight = observation_timings(site, date)

    # specify a time axis
    start_time = sunset - hour
    end_time = sunrise + hour
    delta_t = end_time - start_time
    observe_time = start_time + delta_t * np.linspace(0, 1, 100)

    # get the moon
    moon = get_moon(observe_time)
    illum = moon_illumination(evening_twilight) * 100
    moon_altitude = site.altaz(observe_time, moon).alt.deg

    # plotting
    linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))] * (int(len(targets)/4)+1)
    n = 0
    for obj in targets:
        if min_rank <= obj['rank'] <= max_rank:
            n += 1
            linestyle = linestyles[int((n-1)/n_colors)]
            priority = REV_RANK_INDEX[obj['rank']]
            target_coords = SkyCoord(ra=float(obj['ra'])*u.deg, dec=float(obj['dec'])*u.deg)
            target = FixedTarget(target_coords, name=f'{obj["name"]}({priority[0]})')

            altitude = site.altaz(observe_time, target).alt.deg
            max_alt = np.max(altitude)
            time_max_alt = observe_time[np.argmax(altitude)]

            moon_angle = moon.separation(target_coords)
            angle_at_max = moon_angle[np.argmax(altitude)].deg

            plot_altitude(target, site, observe_time, airmass_yaxis=True, ax=ax,
                          style_kwargs={'linestyle': linestyle, 'linewidth': 1},  # ((7 - obj['rank'])/2) + 0.5},
                          max_altitude=91, min_altitude=15)
            # display moon distance
            text_obj = ax.text(time_max_alt.plot_date, max_alt, f'{angle_at_max:.0f}$^{{\\rm o}}$',
                               ha='center', va='bottom', fontsize='small', zorder=10)
            text_obj.set_bbox({'facecolor': 'white', 'alpha': 0.75, 'linewidth': 0,
                               'boxstyle': 'round, pad=0.0, rounding_size=0.3'})
            # show latest mag
            text_obj = ax.text(time_max_alt.plot_date, max_alt-1, f'{obj["latest mag"]:.1f}',
                               ha='center', va='top', fontsize='small', zorder=10)
            text_obj.set_bbox({'facecolor': 'white', 'alpha': 0.75, 'linewidth': 0,
                               'boxstyle': 'round, pad=0.0, rounding_size=0.3'})

    # lowest viable seeing
    ax.axhline(y=30, color='k', linestyle='--')

    # plot moon track
    ax.plot(observe_time.plot_date, moon_altitude, color='k', linestyle='--', label='Moon')

    # shading
    ax.axvspan(start_time.datetime, sunset.datetime, ymin=0, ymax=1, color='grey', alpha=0.4)
    ax.axvspan(sunset.datetime, evening_twilight.datetime, ymin=0, ymax=1, color='grey', alpha=0.2)
    ax.axvspan(morning_twilight.datetime, sunrise.datetime, ymin=0, ymax=1, color='grey', alpha=0.2)
    ax.axvspan(sunrise.datetime, end_time.datetime, ymin=0, ymax=1, color='grey', alpha=0.4)

    # airmass axis
    airmass_ticks = np.concatenate([np.arange(1, 2.1, 0.1), np.arange(2.2, 3.2, 0.2)])
    altitude_ticks = 90 - np.degrees(np.arccos(1 / airmass_ticks))
    air_ax = fig.get_axes()[-1]
    air_ax.set_yticks(altitude_ticks)
    air_ax.set_yticklabels(np.array([f'{t:.1f}' for t in airmass_ticks]))

    # plot gubbins
    ax.legend(fontsize='small', ncol=7, loc='lower center', frameon=True)
    ax.grid(linestyle='--')
    ax.set(title=f'Highest priority = {REV_RANK_INDEX[min_rank]}    Lowest priority = {REV_RANK_INDEX[max_rank]}  ' +
                 f'  Moon Illumination = {illum:.0f}%')
    plt.tight_layout(pad=0.5)

    return fig


def make_schedule(targets, site, date, min_rank=1, max_rank=6):
    fig = plt.figure()
    gs = gspec.GridSpec(1, 1)
    ax = fig.add_subplot(gs[0, 0])

    n_colors = sum([1 for d in plt.rcParams['axes.prop_cycle']])

    # get timings
    hour = TimeDelta(3600.0 * u.s)
    sec = TimeDelta(1.0 * u.s)
    sunset, sunrise, evening_twilight, morning_twilight = observation_timings(site, date)

    # specify a time axis
    start_time = sunset - hour
    end_time = sunrise + hour
    delta_t = end_time - start_time
    observe_time = start_time + delta_t * np.linspace(0, 1, 75)

    hi_res_time = evening_twilight + (morning_twilight - evening_twilight) * np.linspace(0, 1, 200)

    # plotting
    linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1, 1, 1))] * (int(len(targets)/4)+1)
    n = 0
    for obj in targets:
        if min_rank <= obj['rank'] <= max_rank:
            n += 1
            linestyle = linestyles[int((n-1)/n_colors)]
            priority = REV_RANK_INDEX[obj['rank']]
            target = FixedTarget(SkyCoord(ra=float(obj['ra'])*u.deg, dec=float(obj['dec'])*u.deg))

            altitude = site.altaz(hi_res_time, target).alt.degree
            max_time = hi_res_time[np.argmax(altitude)]

            plot_altitude(target, site, hi_res_time, airmass_yaxis=True, ax=ax,
                          style_kwargs={'linestyle': ':', 'color': 'k', 'alpha': 0.5},
                          max_altitude=91, min_altitude=15)

            if 'ob' in obj.keys():
                ob = obj['ob']
                obs_mask = np.logical_and(max_time - ob['exec_time_s'] * sec <= hi_res_time,
                                          hi_res_time <= max_time)
                if obs_mask.any():
                    ax.plot(hi_res_time[obs_mask].plot_date, altitude[obs_mask],
                            color='w', linestyle='-', linewidth=4)
                    ax.plot(hi_res_time[obs_mask].plot_date, altitude[obs_mask],
                            label=f'{obj["name"]}({priority[0]})', linestyle=linestyle, linewidth=3)
            else:
                ax.plot(hi_res_time.plot_date, altitude,
                        linestyle=linestyle, label=f'{obj["nickname"]}({priority[0]})')

    # lowest viable seeing
    ax.axhline(y=30, color='k', linestyle='--')

    # shading
    ax.axvspan(start_time.datetime, sunset.datetime, ymin=0, ymax=1, color='grey', alpha=0.4)
    ax.axvspan(sunset.datetime, evening_twilight.datetime, ymin=0, ymax=1, color='grey', alpha=0.2)
    ax.axvspan(morning_twilight.datetime, sunrise.datetime, ymin=0, ymax=1, color='grey', alpha=0.2)
    ax.axvspan(sunrise.datetime, end_time.datetime, ymin=0, ymax=1, color='grey', alpha=0.4)

    # airmass axis
    airmass_ticks = np.concatenate([np.arange(1, 2.1, 0.1), np.arange(2.2, 3.2, 0.2)])
    altitude_ticks = 90 - np.degrees(np.arccos(1 / airmass_ticks))
    air_ax = fig.get_axes()[-1]
    air_ax.set_yticks(altitude_ticks)
    air_ax.set_yticklabels(np.array([f'{t:.1f}' for t in airmass_ticks]))

    # plot gubbins
    ax.legend(fontsize='small', ncol=7, loc='lower center', frameon=True)
    ax.grid(linestyle='--')
    ax.set(title=f'Highest priority = {REV_RANK_INDEX[min_rank]}    Lowest priority = {REV_RANK_INDEX[max_rank]}')
    plt.tight_layout(pad=0.5)

    return fig


if __name__ == '__main__':
    # handle login
    if os.path.exists('login.json'):
        with open('login.json', 'r') as f:
            payload = json.load(f)
    else:
        with open('login.json', 'w') as f:
            json.dump({
                'login': 'your_username',
                'password': 'your_password'
            }, f, indent=4)
        print('No login file detected, please enter your details into the new placeholder file.')
        payload = {}

    # import ignore list
    if os.path.exists('ignore_list.txt'):
        ignore_list = list(np.genfromtxt('ignore_list.txt', dtype=str).flatten())
    else:
        np.savetxt(
            'ignore_list.txt', np.array(['# list of targets to ignore in planning e.g. AT2021uuh AT2018cow']), fmt='%s'
        )
        print('Generated placeholder ignore list.')
        ignore_list = []

    # directories
    make_directory('graphs/')
    make_directory('outputs/')

    # acquire transients
    with requests.Session() as s:
        # login to marshall
        s.post('https://www.pessto.org/marshall/', data=payload)
        try:
            response = s.get('https://www.pessto.org/marshall/transients?mwl=pending+observation&format=json')
            classification = response.json()

            response = s.get('https://www.pessto.org/marshall/transients?mwl=following&format=json')
            followup = response.json()
        except json.decoder.JSONDecodeError:
            print('Failed to download targets from marshall! Check your login details.')
            classification = []
            followup = []

    all_transients = []

    # remove each transient without priority or on the ignore list
    classification = [o for o in classification
                      if o['priority'] in RANK_INDEX.keys() and o['name'] not in ignore_list
                      and o['classification date'] is None]
    followup = [o for o in followup
                if o['priority'] in RANK_INDEX.keys() and o['name'] not in ignore_list]

    # label each transient with the observation type and provide a rank and nickname
    for obj in classification:
        obj['obs_type'] = 'classification'
        obj['rank'] = RANK_INDEX[obj['priority']]
        obj['nickname'] = strip_name(obj['name'])

        # choose an OB
        ob_index = np.array([ob['mag_constraint'](float(obj['latest mag'])) for ob in CLASSIFICATION_OBS], dtype=bool)
        this_ob = np.array(CLASSIFICATION_OBS)[ob_index][0]
        obj['ob'] = this_ob

    for obj in followup:
        obj['obs_type'] = 'followup'
        obj['rank'] = RANK_INDEX[obj['priority']]
        obj['nickname'] = strip_name(obj['name'])

    all_transients = classification + followup

    # sort targets based on rank
    all_transients = sorted(all_transients, key=lambda x: x['rank'])
    classification = sorted(classification, key=lambda x: x['rank'])
    followup = sorted(followup, key=lambda x: x['rank'])

    # show the ordered targets
    for obj in all_transients:
        ob = ''
        if 'ob' in obj.keys():
            ob = f'{obj["ob"]["description"]}  {obj["ob"]["exp_time_s"]}s'
        print(obj['name'], obj['priority'])

    # target lists for staralt
    staralt_class = np.array([
        np.array([o['name'], o['ra (sex)'], o['dec (sex)']])
        for o in classification if o['priority'] != 'LOW'
    ])

    staralt_follow = np.array([
        np.array([o['name'], o['ra (sex)'], o['dec (sex)']])
        for o in followup
    ])

    np.savetxt('outputs/followup_list.txt', staralt_follow, fmt='%s')
    np.savetxt('outputs/classification_list.txt', staralt_class, fmt='%s')

    la_silla = Observer.at_site('La Silla Observatory')

    # observability testing
    # plot_wrap(make_schedule(
    #     targets=classification,
    #     site=la_silla,
    #     date=DATE,
    #     max_rank=5
    # ), filename='graphs/test_schedule', show=False)

    # generate altitude plots
    plot_wrap(target_altitudes(
        targets=classification,
        site=la_silla,
        date=DATE,
        min_rank=6,
        max_rank=6,
    ), filename=f'graphs/{DATE}_altitude_class_low')

    plot_wrap(target_altitudes(
        targets=classification,
        site=la_silla,
        date=DATE,
        min_rank=2,
        max_rank=4,
    ), filename=f'graphs/{DATE}_altitude_class_high')

    plot_wrap(target_altitudes(
        targets=followup,
        site=la_silla,
        date=DATE,
        min_rank=1,
        max_rank=5
    ), filename=f'graphs/{DATE}_altitude_follow')

    plot_wrap(target_altitudes(
        targets=all_transients,
        site=la_silla,
        date=DATE
    ), filename=f'graphs/{DATE}_altitude_all')
