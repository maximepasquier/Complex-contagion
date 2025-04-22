import csv

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
    #baseline_props_path = f'LTM_networks/ws/1000/8_props.csv'
    baseline_polarization_path = f'final/LTM/1000/ws/16_top50.csv'
    #eval_file_props_path = f'Evaluation/ws/1000/8_props_baseline.csv'
    eval_file_polarization_path = f'LTM_networks/ws/1000/16/0/0/polarization.csv'

    #print("Props files :")
    #compare_csv(baseline_props_path, eval_file_props_path)
    print("Polarization files :")
    compare_csv(baseline_polarization_path, eval_file_polarization_path)
