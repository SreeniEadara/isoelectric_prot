from Bio import SeqIO
from matplotlib import pyplot

def main():
    print("Sreenivas Eadara, 2024")
    print()
    print("Calculates the isoelectric point and charge at pH")
    print(" of a protein provided in .fasta format.")
    print("If fasta has multiple records, the properties are")
    print(" calculated for the whole complex.")
    try:
        records = SeqIO.parse(input("ENTER FILENAME:\n"), "fasta")
    except:
        print("\nInvalid Path.")
        return

    charged_residues = {
        "NT": [],
        "CT": [],
        "D": 0,
        "E": 0,
        "K": 0,
        "R": 0,
        "H": 0
    }

    for record in records:
        # Add N-termini
        match record.seq[0]:
            case 'A':
                charged_residues["NT"].append(9.7)
            case 'R':
                charged_residues["NT"].append(9.0)
            case 'N':
                charged_residues["NT"].append(8.8)
            case 'D':
                charged_residues["NT"].append(9.8)
            case 'C':
                charged_residues["NT"].append(10.8)
            case 'Q':
                charged_residues["NT"].append(9.1)
            case 'E':
                charged_residues["NT"].append(9.7)
            case 'G':
                charged_residues["NT"].append(9.6)
            case 'H':
                charged_residues["NT"].append(9.2)
            case 'I':
                charged_residues["NT"].append(9.7)
            case 'L':
                charged_residues["NT"].append(9.6)
            case 'K':
                charged_residues["NT"].append(9.0)
            case 'M':
                charged_residues["NT"].append(9.2)
            case 'F':
                charged_residues["NT"].append(9.1)
            case 'P':
                charged_residues["NT"].append(10.6)
            case 'S':
                charged_residues["NT"].append(9.2)
            case 'T':
                charged_residues["NT"].append(10.4)
            case 'W':
                charged_residues["NT"].append(9.4)
            case 'Y':
                charged_residues["NT"].append(9.1)
            case 'V':
                charged_residues["NT"].append(9.6)
        
        # add C-termini
        match record.seq[-1]:
            case 'A':
                charged_residues["CT"].append(2.3)
            case 'R':
                charged_residues["CT"].append(2.2)
            case 'N':
                charged_residues["CT"].append(2.0)
            case 'D':
                charged_residues["CT"].append(2.1)
            case 'C':
                charged_residues["CT"].append(1.7)
            case 'Q':
                charged_residues["CT"].append(2.2)
            case 'E':
                charged_residues["CT"].append(2.1)
            case 'G':
                charged_residues["CT"].append(2.3)
            case 'H':
                charged_residues["CT"].append(1.8)
            case 'I':
                charged_residues["CT"].append(2.4)
            case 'L':
                charged_residues["CT"].append(2.4)
            case 'K':
                charged_residues["CT"].append(2.2)
            case 'M':
                charged_residues["CT"].append(2.3)
            case 'F':
                charged_residues["CT"].append(1.8)
            case 'P':
                charged_residues["CT"].append(2.0)
            case 'S':
                charged_residues["CT"].append(2.2)
            case 'T':
                charged_residues["CT"].append(2.6)
            case 'W':
                charged_residues["CT"].append(2.4)
            case 'Y':
                charged_residues["CT"].append(2.2)
            case 'V':
                charged_residues["CT"].append(2.3)

        for char in record.seq:
            match char:
                case "D":
                    charged_residues["D"]+=1
                case "E":
                    charged_residues["E"]+=1
                case "K":
                    charged_residues["K"]+=1
                case "R":
                    charged_residues["R"]+=1
                case "H":
                    charged_residues["H"]+=1

    # Calculate charges as a function of pH
    pH_range = [ x / 1000 for x in range(0, 14000, 1) ]
    charges_at_pH = []
    for pH in pH_range:
        charge = 0
        for i in range(len(charged_residues["NT"])):
            charge += 1 / (1 + 10**(pH - charged_residues["NT"][i]))
        for i in range(len(charged_residues["CT"])):
            charge -= 1 / (1 + 10**(charged_residues["CT"][i] - pH))
        
        charge -= charged_residues["D"] / (1 + 10**(3.9 - pH))
        charge -= charged_residues["E"] / (1 + 10**(4.3 - pH))
        charge += charged_residues["K"] / (1 + 10**(pH - 10.5))
        charge += charged_residues["R"] / (1 + 10**(pH - 12.5))
        charge += charged_residues["H"] / (1 + 10**(pH - 6.0))

        charges_at_pH.append(charge)
    
    upper_bound = 14
    for i in range(len(charges_at_pH)):
        if charges_at_pH[i] > 0:
            continue
        upper_bound = i / 1000
        break
    lower_bound = 0
    for i in range(len(charges_at_pH) - 1, 0, -1):
        if charges_at_pH[i] < 0:
            continue
        lower_bound = i / 1000
        break

    isoelectric_point = 0.5 * (0.5 * (2 * lower_bound - 1) + 0.5 * (2 * upper_bound + 1))
    print("ISOELECTRIC POINT:")
    print(isoelectric_point)

    pyplot.plot(pH_range, charges_at_pH)
    pyplot.title("Charge vs pH")
    pyplot.xlim(0, 14)
    pyplot.savefig("FIG.png")

    user_input_pH_index = int(float(input("Input Desired pH to get charge:\n")) * 1000)
    print()
    print("CHARGE AT " + str(float(user_input_pH_index) / 1000) + ":")
    print(charges_at_pH[user_input_pH_index])

main()
