#Thomas Hilton Johnson III
#3/15/2019
#Applied Cryptography
#Dr. Sengupta
#DNA Encryption Algorithm
#reference :https://stackoverflow.com/questions/22571259/split-a-string-into-n-equal-parts
#reference:https://stackoverflow.com/questions/32554527/typeerror-list-indices-must-be-integers-or-slices-not-str
#reference:https://stackoverflow.com/questions/10059554/inserting-characters-at-the-start-and-end-of-a-string
#reference: https://docs.python.org/2/library/functions.html#bin
#reference: https://www.pythoncentral.io/cutting-and-slicing-strings-in-python/
#https://stackoverflow.com/questions/7396849/convert-binary-to-ascii-and-vice-versa author jfs
import textwrap

def inputStringDNA():#Function for generating String from user input
    string_DNA = input("Type in the DNA String.\n")
    print(string_DNA)
    return string_DNA

def splitString(strArg):#Function for splitting string into separate parts
    separation_factor = eval(input("Please enter the number of characters to keep per string.\n"))#Input of the number of characters that will be allocated to each partition
    list_of_str = textwrap.wrap(strArg, separation_factor)#Breaks the string into a list of substrings, sized determined by the separation_factor variable
    print(list_of_str)
    return list_of_str

"""
convertDNABinary function converts the amino acids strings to binary digits
as specified by ehat is in nucleicDict argument
"""
def convertDNABinary(strDNA, nucleicDict):
    separation = 1#The separation of each character
    list_of_acids = textwrap.wrap(strDNA, separation)#Breaking each character apart and storing within a list
    binary_list = convertDNAList(list_of_acids, nucleicDict)#Converting the characters to their equivalent in the binary pairings specified in the nucleicDict argument
    binary_str_DNA = ""
    for index in range(binary_list):
        binary_str_DNA = binary_str_DNA + index
    return binary_str_DNA

def convertDNAList(list_arg, nucleicDict):
    for index in range(list_arg):
        if index in nucleicDict:
           index = nucleicDict[index]
        else:
            print("Anomaly present.")
    return list_arg

def addCharacter(listArg):#Can add a character to the beginning of ech partition of the of the substrings that make uyp the larger string
    listSize = len(listArg)#Stores the number of partitions (substrings) of the original string
    for index in range(listSize):#The for loop iterates through the list, selecting each substring
        charInput = input("Input character.\n")#Input a character
        listArg[index] = charInput + listArg[index]#Adds an inputted character at the beginning of the string
        index = index+1# iterates index up
    print(listArg)

def main():
    dnaLetters = {"A":"00", "C": "01", "G": "10", "T": "11"}
    inputtedString = inputStringDNA()
    splitted_Str = splitString(inputtedString)
    addCharacter(splitted_Str)


#main()
print("The binary value of yes.\n")
binary_value = bin(int.from_bytes("yes".encode(), "big"))
print(binary_value)
binary_value = int(binary_value, 2)
print(binary_value)
binary_retrieve = binary_value.to_bytes((binary_value.bit_length()+7) // 8, "big").decode()
print(binary_retrieve)
