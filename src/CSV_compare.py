import csv

'''
Compare two CSV files and print if they are identical or different.

Input :
    - csv1 : path to the first CSV file
    - csv2 : path to the second CSV file
    
Output :
    - Print if the CSV files are identical or different
'''

def read_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        csv_reader = csv.reader(file)
        for row in csv_reader:
            data.append(row)
    return data

def compare_csv(file1_path, file2_path):
    file1_data = read_csv(file1_path)
    file2_data = read_csv(file2_path)

    if file1_data == file2_data:
        print("The CSV files are identical.")
    else:
        print("The CSV files are different.")

if __name__ == "__main__":
    csv1 = f'LTM/ws/1000/16/32/0/polarization.csv'
    csv2 = f'LTM/ws/1000/16/32/0/polarization.csv'
    compare_csv(csv1, csv2)
