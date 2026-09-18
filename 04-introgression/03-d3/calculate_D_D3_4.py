import random
import pandas as pd
import scipy

# Function for the ABBA-BABA statistic.
def calculate_D(sp1_genome, sp2_genome, sp3_genome, sp4_genome):
	
	ABBA_count, BABA_count = 0, 0
	sites = ["A", "T", "G", "C"]
	resolveR = ["A", "G"]
	resolveY = ["C", "T"]
	resolveS = ["G", "C"]
	resolveW = ["A", "T"]
	resolveK = ["G", "T"]
	resolveM = ["A", "C"]

	# To change characters as needed, convert strings to lists.
	sp1_genome_l = list(sp1_genome)
	sp2_genome_l = list(sp2_genome)
	sp3_genome_l = list(sp3_genome)
	sp4_genome_l = list(sp4_genome)

	for i in range(len(sp1_genome_l)): 
		if sp2_genome_l[i] == 'R': # Randomly resolve two-nucleotide ambiguous sites.
			sp2_genome_l[i] = random.choice(resolveR)
		if sp2_genome_l[i] == 'Y':
			sp2_genome_l[i] = random.choice(resolveY)
		if sp2_genome_l[i] == 'S':
			sp2_genome_l[i] = random.choice(resolveS)
		if sp2_genome_l[i] == 'W':
			sp2_genome_l[i] = random.choice(resolveW)
		if sp2_genome_l[i] == 'K':
			sp2_genome_l[i] = random.choice(resolveK)
		if sp2_genome_l[i] == 'M':
			sp2_genome_l[i] = random.choice(resolveM)
		if sp3_genome_l[i] == 'R':
			sp3_genome_l[i] = random.choice(resolveR)
		if sp3_genome_l[i] == 'Y':
			sp3_genome_l[i] = random.choice(resolveY)
		if sp3_genome_l[i] == 'S':
			sp3_genome_l[i] = random.choice(resolveS)
		if sp3_genome_l[i] == 'W':
			sp3_genome_l[i] = random.choice(resolveW)
		if sp3_genome_l[i] == 'K':
			sp3_genome_l[i] = random.choice(resolveK)
		if sp3_genome_l[i] == 'M':
			sp3_genome_l[i] = random.choice(resolveM)
		if sp4_genome_l[i] == 'R': 
			sp4_genome_l[i] = random.choice(resolveR)
		if sp4_genome_l[i] == 'Y':
			sp4_genome_l[i] = random.choice(resolveY)
		if sp4_genome_l[i] == 'S':
			sp4_genome_l[i] = random.choice(resolveS)
		if sp4_genome_l[i] == 'W':
			sp4_genome_l[i] = random.choice(resolveW)
		if sp4_genome_l[i] == 'K':
			sp4_genome_l[i] = random.choice(resolveK)
		if sp4_genome_l[i] == 'M':
			sp4_genome_l[i] = random.choice(resolveM)
		if sp2_genome_l[i] not in sites or sp3_genome_l[i] not in sites or sp4_genome_l[i] not in sites:
			continue # If the site is still unresolved, skip it.

		if sp2_genome_l[i] == sp3_genome_l[i] and sp2_genome_l[i] != sp4_genome_l[i] and sp4_genome_l[i] == sp1_genome_l[i]: # ABBA sites
			ABBA_count += 1
		elif sp2_genome_l[i] == sp4_genome_l[i] and sp3_genome_l[i] != sp4_genome_l[i]	and sp3_genome_l[i] == sp1_genome_l[i]: # BABA sites
			BABA_count += 1

	D = (float(ABBA_count) - float(BABA_count))/(float(ABBA_count) + float(BABA_count)) # D statistic
	
	return D

# Function for D3.
def calculate_D3(sp2_genome, sp3_genome, sp4_genome):
	
	D23, D24, D34 = 0, 0, 0
	sites = ["A", "T", "G", "C"]
	resolveR = ["A", "G"]
	resolveY = ["C", "T"]
	resolveS = ["G", "C"]
	resolveW = ["A", "T"]
	resolveK = ["G", "T"]
	resolveM = ["A", "C"]

	# To change characters as needed, convert strings to lists.
	sp2_genome_l = list(sp2_genome)
	sp3_genome_l = list(sp3_genome)
	sp4_genome_l = list(sp4_genome)

	for i in range(len(sp2_genome_l)): # Count number of pairwise differences.
		if sp2_genome_l[i] == 'R': # Randomly resolve two-nucleotide ambiguous sites.
			sp2_genome_l[i] = random.choice(resolveR)
		if sp2_genome_l[i] == 'Y':
			sp2_genome_l[i] = random.choice(resolveY)
		if sp2_genome_l[i] == 'S':
			sp2_genome_l[i] = random.choice(resolveS)
		if sp2_genome_l[i] == 'W':
			sp2_genome_l[i] = random.choice(resolveW)
		if sp2_genome_l[i] == 'K':
			sp2_genome_l[i] = random.choice(resolveK)
		if sp2_genome_l[i] == 'M':
			sp2_genome_l[i] = random.choice(resolveM)
		if sp3_genome_l[i] == 'R':
			sp3_genome_l[i] = random.choice(resolveR)
		if sp3_genome_l[i] == 'Y':
			sp3_genome_l[i] = random.choice(resolveY)
		if sp3_genome_l[i] == 'S':
			sp3_genome_l[i] = random.choice(resolveS)
		if sp3_genome_l[i] == 'W':
			sp3_genome_l[i] = random.choice(resolveW)
		if sp3_genome_l[i] == 'K':
			sp3_genome_l[i] = random.choice(resolveK)
		if sp3_genome_l[i] == 'M':
			sp3_genome_l[i] = random.choice(resolveM)
		if sp4_genome_l[i] == 'R': 
			sp4_genome_l[i] = random.choice(resolveR)
		if sp4_genome_l[i] == 'Y':
			sp4_genome_l[i] = random.choice(resolveY)
		if sp4_genome_l[i] == 'S':
			sp4_genome_l[i] = random.choice(resolveS)
		if sp4_genome_l[i] == 'W':
			sp4_genome_l[i] = random.choice(resolveW)
		if sp4_genome_l[i] == 'K':
			sp4_genome_l[i] = random.choice(resolveK)
		if sp4_genome_l[i] == 'M':
			sp4_genome_l[i] = random.choice(resolveM)
		if sp2_genome_l[i] not in sites or sp3_genome_l[i] not in sites or sp4_genome_l[i] not in sites:
			continue # If the site is still unresolved, skip it.

		if sp2_genome_l[i] != sp3_genome_l[i]:
			D23 += 1
		if sp2_genome_l[i] != sp4_genome_l[i]:
			D24 += 1
		if sp3_genome_l[i] != sp4_genome_l[i]:
			D34 += 1

	D23, D24, D34 = float(D23)/len(sp2_genome), float(D24)/len(sp2_genome), float(D34)/len(sp2_genome) # divide by genome length

	three_pairwise = [(D24-D23)/(D23+D24), (D23-D34)/(D23+D34), abs(D24-D34)/(D24+D34)] # get the three pairwise divergence comparisons 

	three_pairwise_abs = [abs(value) for value in three_pairwise] # get absolute values of comparisons

	index_min = min(range(len(three_pairwise_abs)), key=three_pairwise_abs.__getitem__) # get index of D3 value

	D3 = three_pairwise[index_min]  # get value of D3

	return D3

# Function to estimate variance of D and D3 using block jackknife
def D_D3_block_bootstrap(sp1_genome, sp2_genome, sp3_genome, sp4_genome, n_replicates):
	
	D_estimates = []
	D3_estimates = []
	
	for i in range(n_replicates): # for each bootstrap replicate

		sp1_bootstrapped, sp2_bootstrapped, sp3_bootstrapped, sp4_bootstrapped = [], [], [], [] #to store bootstrapped genomes

		for j in range(100): # for each resampling window
		
			sampling_index = random.randint(0, 32829) # index to slice, matches input phylip

			#sample windows 
			sp1_bootstrapped.append(sp1_genome[sampling_index:sampling_index+10000])
			sp2_bootstrapped.append(sp2_genome[sampling_index:sampling_index+10000])
			sp3_bootstrapped.append(sp3_genome[sampling_index:sampling_index+10000])
			sp4_bootstrapped.append(sp4_genome[sampling_index:sampling_index+10000])
		
		# Concatenate resampling windows.
		sp1_bootstrapped = ''.join(sp1_bootstrapped)
		sp2_bootstrapped = ''.join(sp2_bootstrapped)
		sp3_bootstrapped = ''.join(sp3_bootstrapped)
		sp4_bootstrapped = ''.join(sp4_bootstrapped)

		# Estimate statistics.
		D_estimates.append(calculate_D(sp1_bootstrapped, sp2_bootstrapped, sp3_bootstrapped, sp4_bootstrapped))
		D3_estimates.append(calculate_D3(sp2_bootstrapped, sp3_bootstrapped, sp4_bootstrapped))

	# Check proportion of positive and negative genomic windows.
	pos_window = 0
	neg_window = 0
	for num in D3_estimates:
		if num > 0:
			pos_window += 1
		if num < 0:
			neg_window += 1
	
	# Estimate means.
	mean_D = sum(D_estimates)/float(len(D_estimates))
	mean_D3 = sum(D3_estimates)/float(len(D3_estimates))	

	# Estimate variance.
	D_stdev = ((sum([(x - mean_D)**2 for x in D_estimates]))/len(D_estimates))**(0.5)
	D3_stdev = ((sum([(x - mean_D3)**2 for x in D3_estimates]))/len(D3_estimates))**(0.5)

	return mean_D, mean_D3, D_stdev, D3_stdev, pos_window, neg_window


# Import a list of trios.
df_trios = pd.read_csv("trios4.csv", sep = '\t')
df_trios.columns = ['spA', 'spB', 'spC']

# Import the phylip file.
df_phylip = pd.read_csv("data.phy")
# Split lines with arbitrary whitespace into separate columns.
df_phylip.columns = ['input']
df_phylip[['idx_num', 'seqs']] = df_phylip['input'].str.split(' ', n = 1, expand = True)

# Perform the following for all trios.
for idx in df_trios.index:
	elements = [] # Defining the four rows to use.
	elements.append('Reference')
	elements.append(df_trios.at[idx, 'spA'])
	elements.append(df_trios.at[idx, 'spB'])
	elements.append(df_trios.at[idx, 'spC'])

	# Get variables with only the relevant rows.
	df_subset1 = df_phylip[df_phylip['idx_num'].str.contains(''.join(elements[0]))].reset_index()
	df_subset2 = df_phylip[df_phylip['idx_num'].str.contains(''.join(elements[1]))].reset_index()
	df_subset3 = df_phylip[df_phylip['idx_num'].str.contains(''.join(elements[2]))].reset_index()
	df_subset4 = df_phylip[df_phylip['idx_num'].str.contains(''.join(elements[3]))].reset_index()

	# Specify the four genomes for the various functions.
	sp1_genome = df_subset1.at[0, 'seqs'].strip()
	sp2_genome = df_subset2.at[0, 'seqs'].strip()
	sp3_genome = df_subset3.at[0, 'seqs'].strip()
	sp4_genome = df_subset4.at[0, 'seqs'].strip()

	# Run functions.
	D_stat = calculate_D(sp1_genome, sp2_genome, sp3_genome, sp4_genome)
	D3_stat = calculate_D3(sp2_genome, sp3_genome, sp4_genome)
	bootstraps = D_D3_block_bootstrap(sp1_genome, sp2_genome, sp3_genome, sp4_genome, 1000)

	# Calculate Z-scores and p-values.
	Zscore_D = (D_stat - bootstraps[0]) / bootstraps[2]
	Zscore_D3 = (D3_stat - bootstraps[1]) / bootstraps[3]
	# pvalue_D = scipy.stats.norm.cdf(Zscore_D)
	# pvalue_D3 = scipy.stats.norm.cdf(Zscore_D3)
	pvalue_D_1tail = scipy.stats.norm.sf(abs(Zscore_D))
	pvalue_D_2tail = scipy.stats.norm.sf(abs(Zscore_D)) * 2
	pvalue_D3_1tail = scipy.stats.norm.sf(abs(Zscore_D3))
	pvalue_D3_2tail = scipy.stats.norm.sf(abs(Zscore_D3)) * 2

	# Print output.
	print(elements[1], elements[2], elements[3], D_stat, D3_stat, bootstraps[0], bootstraps[1], bootstraps[2], bootstraps[3], bootstraps[4], bootstraps[5], Zscore_D, Zscore_D3, pvalue_D_1tail, pvalue_D_2tail, pvalue_D3_1tail, pvalue_D3_2tail)