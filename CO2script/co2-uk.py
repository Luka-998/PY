import requests

BASE_URL = "https://api.carbonintensity.org.uk/intensity" #no available API for serbia in the meantime

def fetch_last_half_hour() -> dict:
    response = requests.get(BASE_URL, timeout=15)
    response.raise_for_status()
    data = response.json()
    return data["data"][0]["intensity"]

def print_intensity_data():
    try:
        intensity_data = fetch_last_half_hour()
        print(f"CO2 Intensity (last half hour): {intensity_data['actual']} gCO2/kWh")
    except requests.RequestException as e:
        print(f"An error occurred while fetching data: {e}")

if __name__ == "__main__":
    print_intensity_data()
