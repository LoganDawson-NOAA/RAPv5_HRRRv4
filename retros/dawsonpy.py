import re, os


# Function to determine which machine we're on
# Code from Mallory Row (NCEP/EMC)
def get_machine():
    hostname = os.environ['HOSTNAME']
    if "machine" in os.environ:
        machine = os.environ['machine']
    else:
        if "MACHINE" in os.environ:
            machine = os.environ['MACHINE']
        else:
            theia_match  = re.match(re.compile(r"^tfe[0-9]{2}$"), hostname)
            hera_match   = re.match(re.compile(r"^hfe[0-9]{2}$"), hostname)
            tide_match   = re.match(re.compile(r"^t[0-9]{2}a[0-9]{1}$"), hostname)
            gyre_match   = re.match(re.compile(r"^g[0-9]{2}a[0-9]{1}$"), hostname)
            surge_match  = re.match(re.compile(r"^slogin[0-9]{1}$"), hostname)
            luna_match   = re.match(re.compile(r"^llogin[0-9]{1}$"), hostname)
            mars_match   = re.match(re.compile(r"^m[0-9]{2}[a-z]{1}[0-9]{1,2}$"), hostname)
            venus_match  = re.match(re.compile(r"^v[0-9]{2}[a-z]{1}[0-9]{1,2}$"), hostname)
            mars_match2  = re.match(re.compile(r"^m[0-9]{2}[a-z]{1}[0-9]{1,2}f$"), hostname)
            venus_match2 = re.match(re.compile(r"^v[0-9]{2}[a-z]{1}[0-9]{1,2}f$"), hostname)
            mars_match3  = re.match(re.compile(r"^m[0-9]{2}[a-z]{1}[0-9]{1,2}.ncep.noaa.gov$"), hostname)
            venus_match3 = re.match(re.compile(r"^v[0-9]{2}[a-z]{1}[0-9]{1,2}.ncep.noaa.gov$"), hostname)
            if hera_match:
                machine = 'HERA'
            elif tide_match or gyre_match:
                machine = 'WCOSS'
            elif luna_match or surge_match:
                machine = 'WCOSS_C'
            elif mars_match or venus_match or mars_match2 or venus_match2 or mars_match3 or venus_march3:
                machine = 'WCOSS_DELL_P3'
            else: 
                print("Cannot find match for "+hostname)
                exit(1)

    return machine, hostname



# Function to submit batch scripts
def submit_job(job_script):
    print('Submitting '+job_script)
    os.system('bsub < '+job_script)



# Function to do multiple replacements using re.sub
# From: https://stackoverflow.com/questions/15175142/how-can-i-do-multiple-substitutions-using-regex-in-python
def multiple_replace(dict, text):
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)

    if __name__ == "__main__":

        text = "Larry Wall is the creator of Perl"

        dict = {
            "Larry Wall" : "Guido van Rossum",
            "creator" : "Benevolent Dictator for Life",
            "Perl" : "Python",
        }

        print multiple_replace(dict, text)



# Function to determine which domains to plot 
def get_domains(case, model_str):
    splains_cases   = ['apr30','may01','may02','may05','may07','may17','may18','may20','may23','may29']
    cplains_cases   = ['may04','may05','may06','may17','may21','may22','may23','may26','may27','may28']
    nplains_cases   = ['may15']
    midwest_cases   = ['may19','may22','may26','may27','may28','may29']
    northeast_cases = ['may19','may23','may26','may28','may29']
    southeast_cases = ['may04','may05','may13']

    if str.upper(model_str) != 'HRRR-AK' and str.upper(model_str) != 'RAP-AK':
        domains = ['conus']
        if str.lower(case) in splains_cases:
            domains.extend(['splains'])
        elif str.lower(case) in nplains_cases:
            domains.extend(['nplains'])
        elif str.lower(case) in midwest_cases:
            domains.extend(['midwest'])
        elif str.lower(case) in northeast_cases:
            domains.extend(['northeast'])
        elif str.lower(case) in southeast_cases:
            domains.extend(['southeast'])
    else:
        domains = ['alaska']
 
    return domains

