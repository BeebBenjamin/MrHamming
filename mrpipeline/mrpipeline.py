#!/usr/bin/python3.6

# import all python3.6 native modules
import operator
import itertools
import sys
import os
import signal
import argparse

from io import StringIO

# import all non-native modules including: distance, biopython and numpy
import distance
import numpy as np

# import various classes and functions required from modules
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align import AlignInfo


def keyboardinterrupthandler(signal):  # noqa
    """In the event of the user interrupting the program using a 'ctrl + c' command
       or 'keyboard interrupt', this function will print the type of interrupt used; 
       restore the terminal cursor to on if it is off; remove the temp_fasta file if
       one is present; and remove a partially complete distances file.

       Args:
                    signal (int): Signal interrupt type number.

       Returns:
                    N/A
    """
    # print a message to the terminal containing the interrupt code, mainly useful to a developer
    print(""""KeyboardInterrupt (ID: {}) has been caught. Cleaning up created files and resetting local terminal 
    settings...""".format(signal))

    # call the setterm function to toggle terminal cursor state
    setterm_cursor()

    # check if temp file is in directory if it is remove it
    if check_file('temp_fasta.fasta'):
        os.remove('temp_fasta.fasta')

    # check if sample file is in directory if it is remove it
    if check_file(sys.argv[2]):
        os.remove(sys.argv[2])

    # exit the program in a clean state
    exit(0)


def setterm_cursor():
    """In order for a progress indicator in the terminal to work properly the cursor
       needs to be switched off during execution of the program.  In the event of an
       error or a keyboard interrupt and or at the end of the program execution the
       state of the cursor needs to be restored so that the user doesn't experience
       any degration in experience in the current terminal session.  If the cursor
       state is currently off or on this function will toggle the state to it's
       opposite.

       Args:
                    N/A

       Returns:
                    N/A
"""
    # if terminal cursor is off set to on
    if os.system('setterm -cursor off'):
        os.system('setterm -cursor on')

    # if terminal cursor is on set to off
    if os.system('setterm -cursor on'):
        os.system('setterm -cursor off')


def check_file(filepath):
    """During program executuion it is necessary to check if a file exists so that there are no errors
       if a file is found to not exist.  This function checks if a file exists and returns True if it 
       does and False if it doesn't. 

       Args:
                    filepath (str): The file path of the file to be checked

       Returns:
                    exists (bool): If the file exists set to 'True' if not set to 'False'
    """
    # define a variable to contain the result of file check
    exists = os.path.isfile(filepath)

    # return a True or a False to the user
    return exists


def hamming(str1, str2):
    """In order to calculate the percentage divergence between a pair of strings
       the hamming distance needs to be calcultated.  This function takes two string 
       objects and compares them returning the hamming distance i.e. the sum of the raw 
       differences between two strings of equal length.

       Args:
                    str1 (str): First DNA sequence in pairwise comparison
                    str2 (str): Second DNA sequence in pairwise comparison

       Returns:
                    no_of_dif (int): Number of differences between strings
    """
    # check strings are the same length, assert is used because it should always be true
    assert len(str1) == len(str2)

    # calculate hamming distance using imap because this is slightly faster than other techniques
    ne = str.__ne__
    ne = operator.ne
    no_of_dif = sum(map(ne, str1, str2))

    # return hamming distance
    return no_of_dif


def count_gaps(sequence):
    """In order to calculate the correct percent divergence between two strings
       where there are gaps present, the gaps need to be counted in order to be deducted
       from the total number of differences.  This function takes a sequence string and 
       returns the number of instances of the gap character ('-') present.

       Args:
                    sequence (str): Sequence string to have gap characters counted

       Returns:
                    no_of_dif (int): Number of differences between strings
    """
    # define gaps as zero so variable can be reused globally
    gaps = 0

    # loop through the sequence character by character and count gaps
    for a in range(0, len(sequence)):
        # if the character is a gap character then increment the gap variable by 1
        if sequence[a] == '-':
            gaps += 1

    # return the number of gaps present
    return gaps


def consensus_seq(str1, str2):
    """In order to calculate the correct percent divergence between two strings
       with gaps present a custom function needs to be used because the functionality 
       required i.e. taking gaps into consideration isn't supplied natively with Biopython
       Therefore this function takes two strings and compares them character by character
       searching for gaps ('-') or abiguous base pairs ('N'). It then builds a consensus string 
       containing the characters that are identical between them, with gaps indicated as a '-' 
       and 'N's where characters differ.

       Note:
                    Currently this function only looks for '-' or 'N' it is possible the user will not
                    provide clean sequences that contain other abiguities.

       Args:
                    str1 (str): First DNA sequence to be used in consensus pair
                    str2 (str): Second DNA sequence to be used in consensus pair

       Returns:
                    consensus_seq (str): The consensus string of the two DNA sequences
    """
    # define an empty list to put the consensus characters into
    consensus_string = []

    # iterate through the length of the two strings with enumeration to use as index values in each string
    for a, b in zip(range(0, len(str1)), range(0, len(str2))):
        # if the characters at the same index between each string are not equal then test for a gap ('-') character
        if str1[a] != str2[b]:
            # noqa if the character at this index for either string is a gap ('-') character then append a gap ('-') to the consensus
            if str1[a] == '-' or str2[b] == '-':
                consensus_string.append('-')
            # noqa otherwise the character at this location must be a different base pair, so append an abiguous base ('N') to the consensus
            else:
                consensus_string.append('N')
        # otherwise the character at this location must be the same character, so append this to the consensus
        else:
            consensus_string.append(str1[a])

        # join the consensus_string (list) into a string variable
        consensus_seq = ''.join(consensus_string)

    # return the consensus sequence
    return consensus_seq


def compare_gaps(str1, str2):
    """In order to calculate the correct percent divergence between two strings
       with gaps present the locations of the gaps in each sequence compared to the
       consensus sequence i.e. a genuine gap needs to be calculated.  This function
       takes two strings and compares if a gap is a missing base pair or a gap across
       a pairwise comparison returning a 'score' matrix in the form of a tuple containing
       the counts of the three different types of gap.  

       Args:
                    str1 (str): First DNA sequence to be used in gap analysis
                    str2 (str): Second DNA sequence to be used in gap analysis

       Returns:
                    gaps_tupple (tupple): A tupple containing the integer counts of the gap types
    """
    # define integer counters on separate lines becuase integers are not iterable
    gaps1 = 0
    gaps2 = 0
    gaps3 = 0

    # iterate through both strings at the same time to make them available for character by character comparisons
    for bp1, bp2 in zip(str1, str2):

        # if base pair in str1 is a gap and is the same across sequences increment gaps1
        if bp1 == "-" and bp1 == bp2:
            gaps1 += 1
        # if base pair in str1 is a gap and is not the same across sequences increment gaps2
        if bp1 == "-" and bp1 != bp2:
            gaps2 += 1
        # if base pair in str2 is a gap and is not the same across sequences increment gaps3
        if bp2 == "-" and bp2 != bp1:
            gaps3 += 1

        # define a tuple using the gap type integer counts
        gaps_tup = [gaps1, gaps2, gaps3]

    # return the gap type tuple
    return gaps_tup


def restart_line():
    """In order to print an updating varible to the command prompt without using a new line
       python needs to write an '\r' to the buffer and print it to the terminal.  This function
       writes '\r' to the stdout buffer and then writes it to the terminal by using 'flush'. 

       Args:
                    N/A

       Returns:
                    N/A
    """
    # write a '\r' to the stdout buffer so the terminal will print to the same line
    sys.stdout.write('\r')
    # flush the stdout buffer so it is printed to the terminal
    sys.stdout.flush()


def mask_sequence(mask, seq):
    """In order to calculate the correct percent divergence between two strings
       with gaps present the locations of the gaps in each sequence compared to the
       alignment (not pairwise) consensus sequence i.e. a genuine gap needs to be calculated.  
       This functiontakes two strings and compares if a gap is a missing base pair or a gap across
       a pairwise comparison returning a 'score' matrix in the form of a tuple containing
       the counts of the three different types of gap.  

       Note:
                    This function is only used if the 'p' flag is present in the arguments passed to the script

       Args:
                    mask (str): Consensus mask that was generated for the whole alignment
                    seq (str): DNA sequence to be masked

       Returns:
                    new_string (str): The DNA sequence masked by the consensus mask

                    or

                    seq (str): The original DNA sequence if it happens to be the same as the mask (unlikely)
    """
    # define a string to contain mutable string (seq is an immutable biopython object)
    new_string = seq
    # if the mask is not equal to the sequence, then iterate through mask
    if mask != seq:
        # iterate through each character in mask using enumerate to create a string index
        for a, char in enumerate(mask):
            # if the character is an ambiguous base pair ('N') or a gap ('-') then replace both with ('-') in sequence
            if char == "N" or char == "-":
                # define a list variable by convert string to list so that the string can be indexed
                new_string = list(new_string)
                # change the character at index to a gap variable
                new_string[a] = '-'

        # join list back into a string object
        new_string = "".join(new_string)

        # return the masked sequence
        return new_string
    else:
        # return the unmodified sequence
        return seq


def main(args) -> bool:
    mode = "normal"

    # print spacer so the output is easier to read
    print("")

    # toggle off the terminal cursor so that the terminal can be updated with an incremental progress counter
    setterm_cursor()

    # check if user has supplied at least 2 arguments because there needs to be an input file and an output file
    if args.outputfile:
        # create a file to output raw distances to
        try:
            raw_distance_file = open(args.outputfile, "w")
        except IOError:
            print(f"\nCannot create file check file path exists!\n")
            return False

    # noqa check sys args for the partial flag to initate a mode shift (used when you want like for like when there is a partial sequence in the mix)
    if args.partial:
        mode = "partial"
        raw_distance_file = open(args.outputfile, "w")

    # grab user defined input file from terminal
    fasta_file = args.fastafile

    # parse the fasta file using biopython to create iterable fasta object
    sequences = SeqIO.parse(open(fasta_file), "fasta")

    # define a list to append sequence data to
    sequence_list = []

    if mode == "partial":
        # align all sequences
        alignment = AlignIO.read(open(fasta_file), "fasta")

        # grab consensus sequence and use as mask
        summary_align = AlignInfo.SummaryInfo(alignment)
        mask = summary_align.gap_consensus(ambiguous='-')

        # define gaps matrix as a 2 dimensional array of zeros
        gaps = np.zeros(shape=(len(alignment), len(mask)))

        # for each sequence in the alignment  where there is not a gap change zero in corresponding gaps matrix to a 1
        for a, fasta in enumerate(alignment):
            for i, base in enumerate(fasta):
                if base == '-':
                    gaps[a][i] = 0
                else:
                    gaps[a][i] = 1

        # define a mask_check variable by summing gaps matrix for each column
        mask_check = gaps.sum(axis=0)

        # noqa define a new_mask variable to store a mask of gaps highlighting the abiguous bases which aren't gaps across sequences
        new_mask = ''

        # noqa compare the mask_check variable the mask, if a score equal to len(alignment) corresponds to a place where there is currently a gap, change value by index in the mask to a ?
        for i, base in enumerate(mask):
            if mask_check[i] == len(alignment) and base == '-':
                new_mask = new_mask + "N"
            else:
                new_mask = new_mask + base

        # replace old mask with new mask set to only mask out consensus gaps
        mask = new_mask

        # define sequences variable and fill by parsing the fasta file
        sequences = SeqIO.parse(open(fasta_file), "fasta")

        # find the number of sequences by counting to variable call no_seqs
        for i, sequence in enumerate(sequences):
            no_seqs = i + 1

        # reparse the fasta file to reset sequences
        sequences = SeqIO.parse(open(fasta_file), "fasta")

        # noqa modify each sequence in sequences applying the consensus gaps mask to not compare differences across gap regions
        for i, record in enumerate(sequences):
            sequence_list.append(record.id + ',' + mask_sequence(mask, record.seq))

            # check if number of sequences is greater than 2, if it is initiate progress feedback
            if no_seqs > 2:
                # define a variable that estimates percentage of sequences masked
                percent_count = "%.2f" % round(i / float(no_seqs - 1) * 100, 2)

                # update the terminal with a progress indicator for how many of the sequences have been masked
                sys.stdout.write("Masking Sequences: " + str(percent_count) + "% complete")
                sys.stdout.flush()
                restart_line()

    # continue and don't apply a consensus gaps mask to the sequences
    else:
        sequences = SeqIO.parse(open(fasta_file), "fasta")
        for record in sequences:
            sequence_list.append(f"{record.id}, {record.seq}")

    # noqa create an array containing all the possible pairwise combinations of all the sequences in the sequence list, check if partial flag used
    pairwise_sequences = itertools.combinations(sequence_list, 2)

    # count the number of pairwise comparisons and put into a variable
    for i, item in enumerate(itertools.combinations(sequence_list, 2)):
        no_pairwise_seqs = i + 1

    # for each pair_wise comparison between sequences
    for i, item in enumerate(pairwise_sequences):

        # create an output file like class that can be written to and read from
        output = StringIO()

        # format the two strings for comparison into fasta format in two string variables
        str1 = f">{item[0].split(',')[0]}\n{item[0].split(',')[1]}\n"
        str2 = f">{item[1].split(',')[0]}\n{item[1].split(',')[1]}\n"

        # create a temporary fasta file on the hard drive
        temp_fasta = open("temp_fasta.fasta", "w")

        # write the fasta formatted sequences to the output class
        output.write(f"{str1}\n")
        output.write(f"{str2}\n")

        # store the value of the output class to a variable
        contents = output.getvalue()

        # close the output class removing it from memory
        output.close()

        # write the contents of the variable to the temporary fasta file
        temp_fasta.write(contents)

        # close the fasta file removing it from memory and saving the changes
        temp_fasta.close()

        # grab temp file name and put into a variable
        temp_fasta = "temp_fasta.fasta"

        # try creating an alignment using the sequences in the temp file
        try:
            alignment = AlignIO.read(open(temp_fasta), "fasta")
        # if unable to align sequences print error message suggesting this error is due to differences in sequence length
        except ValueError:
            print(f"\nCould not align sequences, please check they are all identical in length!\n")
            return False

        # grab the consensus sequence from the pairwise aligment
        consensus = consensus_seq(alignment[0].seq, alignment[1].seq)

        # count the gaps in each sequence including the consensus
        gaps3 = float(count_gaps(alignment[1].seq))
        gaps2 = float(count_gaps(alignment[0].seq))
        gaps1 = float(count_gaps(consensus))

        # use the hammings distance method from the distance module in both directions for both sequences
        dist1 = distance.hamming(alignment[0].seq, alignment[1].seq)

        # compare the locations of gaps between the pairwise seuqneces and their consensus sequence to calculate total genuine gaps
        total_gaps = compare_gaps(alignment[0].seq, alignment[1].seq)

        # calculate the properdistance between sequences excluding gaps
        answer = float(dist1) - (total_gaps[1] + total_gaps[2])

        # check if number of pairwise comparisons is greater than 1, if it is initiate progress feedback
        if no_pairwise_seqs > 1:

            # convert the number of distances into a floating point number that can be converted to a percentage difference later
            percent_count = "%.2f" % round(
                i / float(no_pairwise_seqs - 1) * 100, 2)

            # write an update percentage to the command line to indicate progress through the different combinations
            sys.stdout.write("Distance calculations: " + str(percent_count) + "% complete")
            sys.stdout.flush()
            restart_line()

        # if any gaps in consensus then calculate % based upon length after complete deletion
        if gaps1 > 0:
            # write the pairwise distances to an output file to be used later
            raw_distance_file.write(alignment[0].id + " " + alignment[1].id + " " + str(answer/(len(consensus)-gaps1)) + "\n")
        # calculate % based upon full sequence length
        else:
            # write the pairwise distances to an output file to be used later
            raw_distance_file.write(alignment[0].id + " " + alignment[1].id + " " + str(answer/len(consensus)) + "\n")

    print("")

    # report to the user where they can find their data
    print(f"\nPairwise raw distances written to the following file that can be found in the parent directory: {raw_distance_file.name}\n")

    # close the file containing the pairwise distances
    raw_distance_file.close()

    # remove temp file and clean up
    os.remove('temp_fasta.fasta')

    # toggle on the terminal cursor
    setterm_cursor()


# if script is running from this file run main() and register the keyboard interrupt handler
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculates Hammings distances between species groups in a fasta '
                                                 'alignment')
    parser.add_argument("-ff", "--fastafile", type=str, required=True, help="The alignment you wish to analyse")
    parser.add_argument("-of", "--outputfile", type=str, default="output.dis", help="The path/name of the output file, "
                                                                                    "default is 'output.dis'")
    parser.add_argument("-p", "--partial", action='store_true', help="Ignore gaps and analyse datasets containing "
                                                                     "partial sequences")
    args = parser.parse_args()

    # register keyboard interrupt handler
    signal.signal(signal.SIGINT, keyboardinterrupthandler)

    try:
        success = main(args)
    except BaseException:
        print("Unhandled exception in query builder")
        raise

    if success:
        sys.exit()
    else:
        sys.exit(1)
