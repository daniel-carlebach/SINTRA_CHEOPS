using SINTRA_CHEOPS

# Set FITS and TLE directories to location of existing local data
SINTRA_CHEOPS.fits_dir = "~/cheops/data/fits"
SINTRA_CHEOPS.tle_dir = "~/cheops/data/tle"

# Set email and password for space-track.org
SINTRA_CHEOPS.spacetrack_email = "user@email.com"
SINTRA_CHEOPS.spacetrack_password = "password123"

# Download most recent CHEOPS visits and SatCat
download_satcat()
download_cheops_visits()

# Download images and TLEs for a new file key
file_key = "CH_PR100018_TG045405_V0300"
download_images(file_key)
download_tles(file_key)
