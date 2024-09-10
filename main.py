class Website:
    def __init__(self, url, founding_year, free_to_use):
        self.url = url
        self.founding_year = founding_year
        self.free_to_use = free_to_use

def website_info(website):
    print("URL:", website.url)
    print("Founding year:", website.founding_year)
    print("Free to use:", website.free_to_use)

GitHub = Website("https://github.com", 2008, True)

website_info(GitHub)
