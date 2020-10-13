#!/usr/bin/env python3.6

import sys

[file_name] = sys.argv[1:]

def printHisto(file_name):
	file_handle = open(file_name, "rt")
	numbers = []
		
	for line in file_handle:
		line = line.rstrip('\n')
	
		numbers = list(line.split(","))
		numbers = [int(i) for i in numbers]
		break
	
	length = len(numbers)
	percetage = [0.0]*length
	normal1 =  [0.0]*length
	normal2 = [0.0]*(length-1)
	for i in range(length):
			
		normal1[i] = numbers[i]/(i+1)
		if i < (length-1):
			normal2[i] = numbers[i]/(i+1)
					
	s = sum(numbers)
	s1 = sum(normal1)
	s2 = sum(normal2)
		
	for i in range(length):
		percetage[i] = (numbers[i]/s) * 100
		normal1[i] = (normal1[i]/s1) * 100
		if(i<(length -1)):
			normal2[i] = (normal2[i]/s2) * 100
	print("distribution_of_neighbors:", numbers)
	print("percetage:", percetage)
	print("normal all: ", normal1)
	print("normal -1: ", normal2)
	
printHisto(file_name)

def printHisto(file_name):
	file_handle = open(file_name, "rt")
	numbers = []
		
	for line in file_handle:
		line = line.rstrip('\n')
	
		numbers = list(line.split(","))
		numbers = [int(i) for i in numbers]
		break
	
	length = len(numbers)
	percetage = [0.0]*length
	normal1 =  [0.0]*length
	normal2 = [0.0]*(length-1)
	for i in range(length):
			
		normal1[i] = numbers[i]/(i+1)
		if i < (length-1):
			normal2[i] = numbers[i]/(i+1)
					
	s = sum(numbers)
	s1 = sum(normal1)
	s2 = sum(normal2)
		
	for i in range(length):
		percetage[i] = (numbers[i]/s) * 100
		normal1[i] = (normal1[i]/s1) * 100
		if(i<(length -1)):
			normal2[i] = (normal2[i]/s2) * 100
	print("distribution_of_neighbors:", numbers)
	print("percetage:", percetage)
	print("normal all: ", normal1)
	print("normal -1: ", normal2)