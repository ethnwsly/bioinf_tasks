
with open('spidy_cdsinput.txt', 'r') as input_file:
	read_file = input_file.read()

	new_object = []
	for line in read_file:
		new_object.append(line.strip('\n'))

	print new_object

	input_file.close()


with open('spidy_cdsinput.txt', 'w') as output_file:

	output_file.write('')
	output_file.write('>CDS'+'\n')
	for line in new_object:
		output_file.write(line)
	output_file.close()
