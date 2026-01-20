# Import LeR
from ler import LeR

# Initialize LeR with default settings
ler = LeR()

# Calculate unlensed CBC statistics
unlensed_param = ler.unlensed_cbc_statistics(size=10000, batch_size=10000, resume=False)

# Calculate the detection rate and extract detectable unlensed events
rate_unlensed, unlensed_param_detectable = ler.unlensed_rate()

# Calculate lensed CBC statistics
lensed_param = ler.lensed_cbc_statistics(size=10000, batch_size=10000, resume=False)

# Calculate the detection rate for lensed events
rate_lensed, lensed_param_detectable = ler.lensed_rate()


