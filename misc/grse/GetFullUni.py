## TODO:
# add alliance information for each planet

import time
import requests
from bs4 import BeautifulSoup
import sys
from sys import argv
# User Defined
import random
import pandas as pd

def UUIDcreatePart(length):
    uuidpart = ""
    for _ in range(length):
        uuidchar = format(random.randint(0, 255), 'x')
        if len(uuidchar) == 1:
            uuidchar = "0" + uuidchar
        uuidpart += uuidchar
    return uuidpart

def createUUID():
    return (
        UUIDcreatePart(4) + '-' +
        UUIDcreatePart(2) + '-' +
        UUIDcreatePart(2) + '-' +
        UUIDcreatePart(2) + '-' +
        UUIDcreatePart(6)
    )

# Generate and print a UUID
uuid = createUUID()
print(uuid)
print(len(uuid))




# Function to create a web session
def createWebSession():
    session = requests.Session()
    return session

# Function to handle the login using the web session
def doLogin(session):
    # Display a loading message
    print('Verifying...')

    # Construct the URL for the login page
    site_url = site + 'index.php' + '?page=login'

    # Construct the data payload
    data = {
        'uuid': uuid,
        'login': userName,
        'pass': userPassword
    }

    try:
        # Send a POST request using the session
        response = session.post(site_url, data=data, timeout=10)
        response_text = response.text
        
        # Handle the response from the server
        if response_text == '1':
            loadPage(session, site + 'index.php?page=overview')
            print('Toolbar shown')
        else:
            print('Incorrect data entered!')

    except requests.RequestException as e:
        # Handle request errors
        print('Login Error\n', e)

    # Hide the loading message
    print('Verifying... Done.')

def loadPage(session, page_url):
    try:
        # Send a GET request using the session
        response = session.get(page_url, timeout=10)

    except requests.RequestException as e:
        # Handle request errors
        print('Error while loading the page\n', e)



def parse_short_number(short_number):
    if short_number[-1] == 'K':
        return float(short_number[:-1]) * 1000
    elif short_number[-1] == 'M':
        return float(short_number[:-1]) * 1000000
    else:
        return float(short_number)
    
def get_debris_count(planet_string):
    metal, crystal = 0, 0
    if "Debris field" in planet_string:
        planet_string = planet_string.replace('\xa0', ' ')
        counts = planet_string.split("Debris field")[1].strip().split(" ")
        metal = parse_short_number(counts[0])
        crystal = parse_short_number(counts[2])

    return metal, crystal


### MAIN ###


uuid = "ac47c6f8-3d8a-5588-4ecf-8595a3be0131"
userName = "khan"
userPassword = "Asdfg12345"
site = "https://u6.mmorts.io/"
page_url =  "https://u6.mmorts.io/game.php?page=galaxy&mode=1"
out_csv = argv[1] # output csv file
#previous_csv = argv[2] # previous csv file to compare to

headers = {
        'Accept': 'text/html, */*; q=0.01',
        'User-Agent': 'Mozilla/5.1 (Linux; U; Android; en-us) Android-Phoenix/1.2 AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.71 Safari/547.36',
        'Accept-Encoding': 'gzip, deflate',
        'Accept-Language': 'en-US',
        'X-Requested-With': 'com.munkeez.space'
    }

webSession = createWebSession ()

doLogin(webSession)

# add alliance information for each planet
# 'game.php?page=messages&amp;mode=display&amp;dsp=1&amp;espioopen'
base_url = "https://u6.mmorts.io/game.php?page=statistics"
body = {"uuid":"a156078a-6913-0089-f458-730fa0dfc4c8"
}
response = webSession.get(base_url, headers=headers, data=body)
soup = BeautifulSoup(response.content, 'html.parser')
# Find the table element
tables = soup.find_all('table')
table = tables[1]
# Find the tbody element
tbody = table.find('tbody')
player_dict = {}
# Iterate through rows in tbody
for row in tbody.find_all('tr'):
    cells = row.find_all('td')

    # Extract data from cells
    rank = cells[0].text.strip()
    plus_minus = cells[1].text.strip()  # This might include "+" or "-"
    player = cells[2].text.strip()
    status = cells[3].text.strip()
    alliance = cells[4].text.strip()
    points = cells[5].text.strip()

    # Print or process the extracted data as needed
    player_dict[player] = {'rank': rank, 'alliance': alliance, 'total_points': points}

print(player_dict)
### Scanning the universe ###

galTemp = []
for galaxyNumber in [1, 2, 3]:
    for systemNumber in range(1, 500):
        remoteFileLocation = "https://u6.mmorts.io/game.php?page=galaxy&mode=1"
        print(f"Inspecting {galaxyNumber}:{systemNumber}")
        body = {"galaxy":galaxyNumber,
                "system":systemNumber,
                "uuid":"a156078a-6913-0089-f458-730fa0dfc4c8"}
        response = webSession.post(remoteFileLocation, headers=headers, data=body)
        soup = BeautifulSoup(response.content, 'html.parser')
        planet_info = soup.select('.g_right')
        systemTemp = []

        for incr, planet in enumerate(planet_info, start=1):
            coords = f"{galaxyNumber}:{systemNumber}:{incr}"
            fullString = planet.get_text(strip=True).split('\n')[0]
            #print(fullString)
            if fullString == "":
                continue # no planet
            else:
                name = fullString.split('(ranked')[0]
                #print(name)
                inactive = '(i)' in fullString or '(I)' in fullString
                rank = fullString.split()[1].split(")")[0]
                metal_df, crystal_df = get_debris_count(fullString)
                tempVar = {
                    'coords': coords,
                    'galaxy': galaxyNumber,
                    'system': systemNumber,
                    'planet': incr,
                    'name': name,
                    'inactive': inactive,
                    'rank': rank,
                    'metal_df': metal_df,
                    'crystal_df': crystal_df
                }
                systemTemp.append(tempVar)
        time.sleep(round(random.random()*2, 3)) # Sleep between 0 and 4 seconds

        galTemp.extend(systemTemp)


# convert 'inactive' column to 'status', if inactive then 1 else 0 
for planet_info in galTemp:
    planet_info['status'] = "Inactive" if planet_info['inactive'] else ""
    del planet_info['inactive']

# print the collected data in a dataframe
df = pd.DataFrame(galTemp)

# add alliance column
alliances = []
for player in df['name']:
    try:
        alliances.append(player_dict[player]['alliance'])
    except KeyError:
        alliances.append("")
df['alliance'] = alliances
# sort df by crystal_df
df.sort_values(by=['crystal_df'], inplace=True, ascending=False)
df.to_csv(out_csv, index=False)

'''
# compare 2 datarames:
prev_df = pd.read_csv(previous_csv)

# create unique ID for each ow by combining coords and name
df['unique_id'] = df['coords'] + df['name']
prev_df['unique_id'] = prev_df['coords'] + prev_df['name']
# check if there are differences between sets of unique IDs between 2 dataframes
print("Unique planets: Current")
print(set(df['unique_id']) - set(prev_df['unique_id']))
print("Unique planets: Previous")
print(set(prev_df['unique_id']) - set(df['unique_id']))
'''

