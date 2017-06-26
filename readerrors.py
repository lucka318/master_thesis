import numpy as np
import pysam
import re
from collections import defaultdict

class State:
	on_on = 0
	on_off = 1
	off_on = 2
	off_off = 3

def make_freq_csv(freq_dict, read_max, ref_max, out):
	#print(freq_dict)
	freq_cnt = np.zeros(shape=(ref_max + 1, read_max + 1), dtype=float)
	total = 0
	for key, val in freq_dict.items():
		ref_len = key[0]
		read_len = key[1]
		letter = key[2]
		freq_cnt[ref_len][read_len] += val
		total += val
	print(total)
	freq_cnt = np.divide(freq_cnt, float(total))
	with open(out,'w') as f:
		np.savetxt(f, freq_cnt.astype(float), fmt='%2.5f', delimiter=',')

def make_test_csv(regions_dict, read_max, out):
	size = len(regions_dict)
	freq_cnt = np.zeros(shape=(size, read_max + 4), dtype=int)
	cnt = 0
	for key, val in regions_dict.items():
		if(cnt >= size):
			break
		length = key[1] - key[0]
		for i in val:
			freq_cnt[cnt][i] += 1
		freq_cnt[cnt][read_max + 1] = length
		freq_cnt[cnt][read_max + 2] = key[0]
		freq_cnt[cnt][read_max + 3] = key[1]
		cnt += 1
	with open(out,'w') as f:
		np.savetxt(f, freq_cnt.astype(int), fmt='%i', delimiter=',')

def genome_preprocessing(reference_file):

	with open(reference_file) as f:
		content = f.readlines()
	content = [x.strip() for x in content]
	genome=''.join(content[1:])
	return genome

def get_cigar_string(cigar_string):

	cigar_components = re.findall(r'(\d+)([A-Z]{1})', cigar_string)
	cigar = ''
	for comp in cigar_components:
		number = comp[0]
		letter = comp[1]
		for i in range(int(number)):
			cigar+=letter
	return cigar


def readerrors(ref, reads):

	genome = genome_preprocessing(ref) + 'X'
	content = pysam.AlignmentFile(reads)

	freq_dict = {}
	regions_dict = defaultdict(list)

	#track max homopolymer size
	ref_max = 0
	read_max = 0
	for r in content.fetch(until_eof=True):

# 		#ako se read nije namapirao na referencu
 		if r.cigarstring is None:
 			continue

 		cigar = get_cigar_string(r.cigarstring)

		#set up necessary variables
		#current pointers, they show us which positions we are looking at
		reference_pointer = r.pos
		read_pointer = 0
		cigar_pointer = 0

		#length of homopolymers
		read_homopolymer = 0
		ref_homopolymer = 0

		#coordinates of homopolymers
		ref_begin = 0
		ref_end = 0

		#homopolymer letters. if they exist, they should never be different (except if one of them is empty)
		read_letter = ''
		ref_letter = ''

		#read sequence
		read = r.seq + 'X'

		#starting state
		state = State.off_off


		while(cigar_pointer < len(cigar)):

			if(cigar[cigar_pointer] == 'M'):

			#not a read nor a reference has a detected homopolymer
				if(state == State.off_off):
					if(read[read_pointer] == read[read_pointer + 1]):
						read_letter = read[read_pointer]
						read_homopolymer += 1
						state = State.on_off

						if(genome[reference_pointer] == read_letter):
							ref_homopolymer += 1
							ref_letter = read_letter
							state = State.on_on
							ref_begin = reference_pointer

					elif(genome[reference_pointer] == genome[reference_pointer + 1]):
						ref_homopolymer += 1
						ref_letter = genome[reference_pointer]
						state = State.off_on
						ref_begin = reference_pointer

						if(read[read_pointer] == ref_letter):
							read_letter = ref_letter
							read_homopolymer += 1
							state = State.on_on

					reference_pointer += 1
					read_pointer += 1
					cigar_pointer += 1

				# we have a homopolymer in the read, but not in the reference(or ref homopolymer ended)
				elif(state == State.on_off):

					#if homopolymer in read continues
					if(read[read_pointer] == read_letter):
						read_homopolymer += 1
						#if we didn't found a homopolymer in the reference by now, check
						if(ref_homopolymer == 0 and genome[reference_pointer] == read_letter):
							ref_letter = read_letter
							ref_homopolymer += 1
							state = State.on_on
							ref_begin = reference_pointer

						reference_pointer += 1
						read_pointer += 1
						cigar_pointer += 1

					else:
						#here we don't update pointers so we can check them in off_off state
						state = State.off_off

						data = (ref_homopolymer, read_homopolymer, read_letter)
						if(freq_dict.has_key(data)):
							freq_dict[data] += 1
						else:
							freq_dict[data] = 1

						ref_coordinates = (ref_begin, ref_end)
						regions_dict[ref_coordinates].append(read_homopolymer)


						ref_max = max(ref_max, ref_homopolymer)
						read_max = max(read_max, read_homopolymer)

						ref_begin = 0
						ref_end = 0
						read_homopolymer = 0
						ref_homopolymer = 0
						read_letter = ''
						ref_letter = ''

				# we have a homopolymer in the reference, but not in the read (or read homopolymer ended)
				elif(state == State.off_on):

					#if homopolymer in reference continues
					if(genome[reference_pointer] == ref_letter):
						ref_homopolymer += 1
						if(read_homopolymer == 0 and read[read_pointer] == ref_letter):
							read_letter = ref_letter
							read_homopolymer += 1
							state = State.on_on

						reference_pointer += 1
						read_pointer += 1
						cigar_pointer += 1

					else:
						#here we don't update pointers so we can check them in off_off state
						state = State.off_off
						ref_end = reference_pointer

						data = (ref_homopolymer, read_homopolymer, ref_letter)
						if(freq_dict.has_key(data)):
							freq_dict[data] += 1
						else:
							freq_dict[data] = 1

						ref_coordinates = (ref_begin, ref_end)
						regions_dict[ref_coordinates].append(read_homopolymer)

						ref_max = max(ref_max, ref_homopolymer)
						read_max = max(read_max, read_homopolymer)

						ref_begin = 0
						ref_end = 0
						read_homopolymer = 0
						ref_homopolymer = 0
						read_letter = ''
						ref_letter = ''

				elif(state == State.on_on):
					#print("ON_ON")
					# if read homopolymer continues
					if(read[read_pointer] == read_letter and genome[reference_pointer] == ref_letter):
						read_homopolymer += 1
						ref_homopolymer += 1

						reference_pointer += 1
						read_pointer += 1
						cigar_pointer += 1

					elif(read[read_pointer] == read_letter):
						read_homopolymer += 1
						state = State.on_off

						ref_end = reference_pointer

						reference_pointer += 1
						read_pointer += 1
						cigar_pointer += 1

					elif(genome[reference_pointer] == ref_letter):
						ref_homopolymer += 1
						state = State.off_on

						reference_pointer += 1
						read_pointer += 1
						cigar_pointer += 1

					else:
						state = State.off_off

						ref_end = reference_pointer

						data = (ref_homopolymer, read_homopolymer, read_letter)
						if(freq_dict.has_key(data)):
							freq_dict[data] += 1
						else:
							freq_dict[data] = 1

						ref_coordinates = (ref_begin, ref_end)
						regions_dict[ref_coordinates].append(read_homopolymer)

						ref_max = max(ref_max, ref_homopolymer)
						read_max = max(read_max, read_homopolymer)

						ref_begin = 0
						ref_end = 0

						read_homopolymer = 0
						ref_homopolymer = 0
						read_letter = ''
						ref_letter = ''

			elif(cigar[cigar_pointer] == 'I'):

				if(state == State.off_off):
					#print("OFF_OFF")
					if(read[read_pointer] == read[read_pointer + 1]):
						read_letter = read[read_pointer]
						read_homopolymer += 1
						state = State.on_off

					read_pointer += 1
					cigar_pointer += 1
					# we have a homopolymer in the read, but not in the reference(or ref homopolymer ended)
				elif(state == State.on_off):
					#print("ON_OFF")
					if(read[read_pointer] == read_letter):
						read_homopolymer += 1
						read_pointer += 1
						cigar_pointer += 1
					else:
						state = State.off_off

						data = (ref_homopolymer, read_homopolymer, read_letter)
						if(freq_dict.has_key(data)):
							freq_dict[data] += 1
						else:
							freq_dict[data] = 1

						ref_coordinates = (ref_begin, ref_end)
						regions_dict[ref_coordinates].append(read_homopolymer)

						ref_max = max(ref_max, ref_homopolymer)
						read_max = max(read_max, read_homopolymer)

						ref_begin = 0
						ref_end = 0

						read_homopolymer = 0
						ref_homopolymer = 0
						read_letter = ''
						ref_letter = ''

				elif(state == State.off_on):
					#print("OFF_ON")
					if(read_homopolymer == 0 and read[read_pointer] == ref_letter):
						read_letter = ref_letter
						read_homopolymer += 1
						state = State.on_on
					read_pointer += 1
					cigar_pointer += 1

				elif(state == State.on_on):
					#print("ON_ON")
					if(read[read_pointer] == read_letter):
						read_homopolymer += 1
					else:
						state = State.off_on
					read_pointer += 1
					cigar_pointer += 1


			elif(cigar[cigar_pointer] == 'D'):
				if(state == State.off_off):
					#print("OFF_OFF")
					if(genome[reference_pointer] == genome[reference_pointer + 1]):
						ref_homopolymer += 1
						ref_letter = genome[reference_pointer]
						ref_begin = reference_pointer
						state = State.off_on
					reference_pointer += 1
					cigar_pointer += 1
					# we have a homopolymer in the read, but not in the reference(or ref homopolymer ended)
				elif(state == State.on_off):
					#print("ON_OFF")
					if(ref_homopolymer == 0 and genome[reference_pointer] == read_letter):
						ref_letter = read_letter
						ref_homopolymer += 1
						state = State.on_on
						ref_begin = reference_pointer
					reference_pointer += 1
					cigar_pointer += 1

				elif(state == State.off_on):
					#rint("OFF_ON")
					if(genome[reference_pointer] == ref_letter):
						ref_homopolymer += 1
						reference_pointer += 1
						cigar_pointer += 1
					else:
						state = State.off_off

						ref_end = reference_pointer

						data = (ref_homopolymer, read_homopolymer, ref_letter)
						if(freq_dict.has_key(data)):
							freq_dict[data] += 1
						else:
							freq_dict[data] = 1

						ref_coordinates = (ref_begin, ref_end)
						regions_dict[ref_coordinates].append(read_homopolymer)

						ref_max = max(ref_max, ref_homopolymer)
						read_max = max(read_max, read_homopolymer)

						ref_begin = 0
						ref_end = 0

						read_homopolymer = 0
						ref_homopolymer = 0
						read_letter = ''
						ref_letter = ''

				elif(state == State.on_on):
					#print("ON_ON")
					if(genome[reference_pointer] == ref_letter):
						ref_homopolymer += 1
					else:
						state = State.on_off
						ref_end = reference_pointer
					cigar_pointer += 1
					reference_pointer += 1
			else:
				# S
				read_pointer += 1
				cigar_pointer += 1

	return read_max, regions_dict