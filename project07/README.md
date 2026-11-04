# Introduction
A sequencer produces hundreds of millions of short DNA reads, and each one has to be located in a reference genome of billions of bases. Scanning the genome once per read is too slow, so aligners build an **index** of the reference once and then search it quickly for every read. The most widely used short-read aligners, such as Bowtie (Langmead et al., 2009), BWA (Li & Durbin, 2009) and Bowtie 2 (Langmead & Salzberg, 2012), use an FM-index built on the **Burrows–Wheeler Transform (BWT)** (Burrows & Wheeler, 1994; Ferragina & Manzini, 2005).

The BWT of a text is made by adding an end-of-string character `$`, sorting every rotation of the text, and reading off the last column. The transform is a permutation of the text, it groups identical characters into runs, and it is reversible. The property that makes searching (and reversing) possible is the **last-to-first (LF) mapping**: the k-th occurrence of a character in the last column is the same text character as its k-th occurrence in the first column.

Our program does this in three stages:
1. **Build the BWT:** add the end-of-string character `$`, then build the transform two ways: by sorting all rotations (naive), and from a suffix array (efficient). Both give the same string.
2. **Pattern matching with the LF mapping:** build the `count` and `occur` tables, which together give the LF mapping, then run backward search: read the pattern right to left and narrow a range of rows. The suffix array turns the final rows into positions in the text.
3. **Compress it:** run-length encode the BWT string, and decode it back.

# Documentation and Results

### Requirements

-   Python 3.14
-   Jupyter Notebook
-   Python standard library only: `random` (test data) and `collections.Counter`


### Usage

Open `project07.ipynb` in Jupyter and run all cells in order, top to bottom (Kernel → Restart & Run All). Each function's cell ends with its own tests, printed as `result | expected: ...`.\
The notebook is split into:\
Test data: `generate_dna()` (a reproducible random DNA sequence)\
Stage 1 functions: `BWT()`, `suffix_array()`, `BWT_from_suffix_array()`\
Stage 2 functions: `cal_count()`, `cal_occur()`, `update_range()`, `find_match()`\
Stage 3 functions: `run_length_encode()`, `run_length_decode()`\
Demo: `main()` runs pattern search and compression on the DNA sequence

### Project structure

```
project07/
├── project07.ipynb
├── README.md
```

### Data

There is no data file. We tested with:
- the examples in the starter notebook's docstrings (e.g. `banana`, `googol`, `ATGACG`, `AG$CG`)
- a reproducible 150-bp random DNA sequence, `generate_dna(150, seed=42)`, so the outputs are the same every run

### Checks

Pattern-matching tests (the assignment asks for a pattern that occurs once, several times, overlaps itself, and does not occur):

| Case | Query | Result |
|---|---|---|
| Occurs once | `dna[10:20]` = `'AGAAGGAGTA'` | `[10]` |
| Occurs several times | `'GAC'` | `[0, 57, 105, 123, 147]` |
| Occurs several times | `'CATT'` | `[69, 125, 130, 135]` |
| Overlaps itself | `'ana'` in `banana` | `[1, 3]`: the match at 1 (positions 1–3) and the match at 3 (positions 3–5) share position 3 |
| Does not occur | `'TTTT'` | `[]` |
| Character not in the reference | `'NNN'`, and `'xyz'` in `banana` | `[]` |
| Docstring example | `'na'` in `banana` | `[2, 4]` |

Checking the inverse: the notebook has no separate inverse function. Reversibility follows from the LF mapping (see the pseudocode's *Checking the inverse* section). Traced by hand for `banana`, rows 0 → 1 → 5 → 2 → 6 → 3 → 4 rebuild `banana$`.

### Design questions

**Inputs, outputs, and the BWT property each function relies on:**

| # | Function | Inputs → Output | Property it relies on |
|---|---|---|---|
| 1 | `BWT` | text → BWT string | The definition: sort the rotations, take the last column|
| 2 | `suffix_array` | text → sorted suffix positions | With a unique, smallest `$`, suffixes sort in the same order as rotations |
| 3 | `BWT_from_suffix_array` | text, suffix array → BWT string | The last character of each row comes just before its first character in the text, so BWT[i] = text[SA[i] − 1] |
| 4 | `cal_count` | text without `$` → `{char: block start}` | The first column is the text sorted, so each character's rows form one block |
| 5 | `cal_occur` | BWT → `{char: running tally}` | LF needs "which occurrence" a character is (1st, 2nd, …) |
| 6 | `update_range` | range, `count`, `occur`, char → new range | LF mapping applied to the top and bottom of a block of rows |
| 7 | `find_match` | query, reference without `$` → positions | Every match is a prefix of some suffix, and those rows form one block; LF walks backward, so the query is read right to left |
| 8 | `run_length_encode` | BWT → encoded string | The BWT groups identical characters into runs |
| 9 | `run_length_decode` | encoded string → BWT | None of its own (plain run-length decoding); it gives back the exact BWT |

# Pseudocode
The `# Check` lines are the notebook's tests: the docstring examples first, then the reproducible 150-bp random DNA sequence (`dna_seq = generate_dna(150, seed=42)`).

## Setup: test data
0. Generate a random DNA sequence
```
def generate_dna(length = 150, seed = 42):
    #1. Fix the random seed
    # same seed => same sequence every run, so test outputs are reproducible
    set the random seed to seed

    #2. Pick length random bases from A, C, G, T and join them
    return a string of length random choices from "ACGT"

    # Check: generate_dna() -> 150 bp, starting 'GACAGGTACAAGAAGGAGTA...'
```

## Stage 1: Build the BWT
1. BWT (sort rotations)
```
def BWT(string):
    #1. Add the stop sign $ only if it's missing
    # first check: '$' may appear at most once, and only at the end (raise an error if not)
    if string has more than one '$', or a '$' that is not at the end:
        raise an error
    if string does not end with '$':
        string = string + '$'

    #2. Build every rotation: add the current one, then move its first letter to the back
    rotations = empty list
    current = string
    for each character in string:
        append current to rotations
        current = current[1:] + current[0]

    #3. Sort rotations lexicographically
    # .sort() compares letter by letter using ord(): $ = 36, a = 97, b = 98, n = 110
    sort rotations

    #4. Take the LAST character of each sorted rotation
    bwt = join the last character of every rotation
    return bwt

    # Check: BWT('googol')  -> 'lo$oogg'
    #        BWT('banana$') -> 'annb$aa'
    #        BWT(dna_seq)   -> length 151 (= 150 + 1 for '$')
```
2. Suffix array
```
def suffix_array(string):
    #1. Add the stop sign $ only if it's missing
    # first check: '$' may appear at most once, and only at the end (raise an error if not)
    if string has more than one '$', or a '$' that is not at the end:
        raise an error
    if string does not end with '$':
        string = string + '$'

    #2. Make one (suffix, position) pair per starting position
    pairs = empty list
    for i from 0 to len(string) - 1:
        append (string[i:], i) to pairs

    #3. Sort by the suffix (the first item of each pair)
    sort pairs

    #4. Keep only the positions, in sorted order
    positions = the second item of every pair
    return positions

    # Check: suffix_array('googol')  -> [6, 3, 0, 5, 2, 4, 1]
    #        suffix_array('banana$') -> [6, 5, 3, 1, 0, 4, 2]
    #        suffix_array(dna_seq)[0] == 150   (the '$' suffix always sorts first)
```
3. BWT from the suffix array
```
def BWT_from_suffix_array(text, suffix_positions):
    # suffix_positions comes from suffix_array(text)

    #1. Add the stop sign $ only if it's missing (same as BWT and suffix_array)
    # first check: '$' may appear at most once, and only at the end (raise an error if not)
    if text has more than one '$', or a '$' that is not at the end:
        raise an error
    if text does not end with '$':
        text = text + '$'
    bwt = ""

    #2. For each position p in the suffix array, take the character right before p
    # BWT[i] = text[SA[i] - 1]; p = 0 has nothing before it, so wrap to the last char ($)
    for each p in suffix_positions:
        if p == 0:
            char = last character of text
        else:
            char = text[p - 1]
        bwt = bwt + char

    #3. Return the finished BWT
    return bwt

    # Check: BWT_from_suffix_array('banana$', [6, 5, 3, 1, 0, 4, 2]) -> 'annb$aa'
    #        BWT_from_suffix_array('googol$', [6, 3, 0, 5, 2, 4, 1]) -> 'lo$oogg'
    #        BWT_from_suffix_array(dna_seq + '$', suffix_array(dna_seq)) == BWT(dna_seq) -> True
```

## Stage 2: Pattern matching (with the LF mapping)
```
   #LF mapping: the k-th c in the last column is the k-th c in the first column
   #data structures: count[c]      = where c's block starts in the first column  (function 4)
   #                 occur[c][row] = which c it is (1st, 2nd, ...), counting THROUGH row  (function 5)
   #so the k-th c sits in row count[c] + k, and LF(row) = count[c] + occur[c][row]
```
4. Count dictionary
```
def cal_count(string):
    #1. Count how many times each character shows up
    # use Counter
    char_counts = Counter(string)

    #2. Walk through the characters in lexicographical order
    # store the running total FIRST (only smaller characters counted so far),
    # then add this character's own count for the next one
    running_total = 0
    count = empty dictionary
    for each char in sorted(char_counts):
        count[char] = running_total
        running_total = running_total + char_counts[char]

    #3. Return the finished dictionary
    return count

    # Check: cal_count('ATGACG') -> {'A': 0, 'C': 2, 'G': 3, 'T': 5}
    #        cal_count('banana') -> {'a': 0, 'b': 3, 'n': 4}
    #        cal_count(dna_seq)  -> {'A': 0, 'C': 41, 'G': 74, 'T': 114}
    # Note: '$' is not skipped, so call it on the reference WITHOUT '$'
    # this is the "where c's block starts" half of the LF mapping
```
5. Occur dictionary
```
def cal_occur(bwt_string):
    #1. Give each unique character its own list of zeros
    alphabet = sorted unique characters of bwt_string
    occur = for each char in alphabet: a list of len(bwt_string) zeros

    #2. Walk the BWT one position at a time
    # copy every tally from the previous position (position 0 has nothing to copy),
    # then add 1 for the character sitting at this position
    # occur[c][i] = how many c's from position 0 THROUGH position i
    for i from 0 to len(bwt_string) - 1:
        if i > 0:
            for each char in alphabet:
                occur[char][i] = occur[char][i - 1]
        occur[bwt_string[i]][i] = occur[bwt_string[i]][i] + 1

    #3. Return the finished table
    # this is the "which occurrence" half of the LF mapping
    return occur

    # Check: cal_occur('AG$CG')   -> {'$': [0, 0, 1, 1, 1], 'A': [1, 1, 1, 1, 1],
    #                                 'C': [0, 0, 0, 1, 1], 'G': [0, 1, 1, 1, 2]}
    #        cal_occur('annb$aa') -> as in the example above
    #                                (the notebook's docstring shows a different answer)
    #        the last tally of each letter in cal_occur(BWT(dna_seq)) equals its count in the BWT:
    #        {'$': 1, 'A': 41, 'C': 33, 'G': 40, 'T': 36}
```
6. Update the range
```
def update_range(lower, upper, count, occur, a):
    #1. Tally of a just BEFORE the range starts
    # if lower is 0 there is nothing before it, so use 0 (not occur[a][-1])
    if lower == 0:
        tally_before = 0
    else:
        tally_before = occur[a][lower - 1]

    #2. Tally of a at the END of the range
    tally_end = occur[a][upper]

    #3. New lower: jump to a's block in the first column (+1 for the '$' row), skip earlier a's
    # = row of the FIRST a in the range (LF of the top)
    new_lower = count[a] + 1 + tally_before

    #4. New upper: the last a inside the range
    # = row of the LAST a in the range (LF of the bottom)
    new_upper = count[a] + tally_end

    #5. Return the new range
    return (new_lower, new_upper)

    # Check (banana, 'ana' right to left): update_range(0, 6, ..., 'a') -> (1, 3)
    #                                      update_range(1, 3, ..., 'n') -> (5, 6)
    #                                      update_range(5, 6, ..., 'a') -> (2, 3)
    #        DNA, whole table (rows 0–150), one base each: block size = that base's count
    #        'A' -> rows 1–41 (41), 'C' -> rows 42–74 (33), 'G' -> rows 75–114 (40), 'T' -> rows 115–150 (36)
```
7. Find matches
```
def find_match(query, reference):
    #1. Add $, then build the suffix array and the BWT
    text = reference + '$'
    sa = suffix_array(text)
    bwt = BWT_from_suffix_array(text, sa)

    #2. Count from the plain reference, occur from the BWT
    count = cal_count(reference)
    occur = cal_occur(bwt)

    #3. Start with every row as a possible match
    lower = 0
    upper = len(bwt) - 1

    #4. Read the query from right to left, shrinking the range
    # stop if a letter isn't in the reference, or the range becomes empty (lower > upper)
    for each letter in reversed(query):
        if letter is not in count:
            return empty list
        (lower, upper) = update_range(lower, upper, count, occur, letter)
        if lower > upper:
            return empty list

    #5. Turn the remaining rows into positions using the suffix array, then sort
    positions = empty list
    for each row from lower to upper:
        append sa[row] to positions
    return sorted(positions)

    # Check: find_match('ana', 'banana') -> [1, 3]                       (overlaps itself: both matches use position 3)
    #        find_match('xyz', 'banana') -> []
    #        find_match('na', 'banana')  -> [2, 4]
    #        DNA: dna_seq[10:20] = 'AGAAGGAGTA' -> [10]                        (occurs once)
    #             'GAC'    -> [0, 57, 105, 123, 147]                          (several times)
    #             'CATT'   -> [69, 125, 130, 135]                              (several times)
    #             'TTTT'   -> []                                               (does not occur)
    #             'NNN'    -> []                                               (letter not in the reference)
```

## Stage 3: Compress
```
   #format: the count goes BEFORE the character, and a count of 1 is left out
```
8. Run-length encode
```
def run_length_encode(bwt_string):
    #1. Nothing to encode
    if bwt_string is empty:
        return ""

    #2. Start tracking the first character
    encoded = ""
    current_char = first character of bwt_string
    current_count = 1

    #3. Go through the rest of the string
    # same letter => the run keeps going
    # different letter => the run ended: write it out (count first, only if > 1), then track the new letter
    for each char in bwt_string, starting from the second:
        if char == current_char:
            current_count = current_count + 1
        else:
            if current_count > 1:
                encoded = encoded + current_count + current_char
            else:
                encoded = encoded + current_char
            current_char = char
            current_count = 1

    #4. Write out the last run (it never reaches the else above)
    if current_count > 1:
        encoded = encoded + current_count + current_char
    else:
        encoded = encoded + current_char

    return encoded

    # Check: run_length_encode('aaabbc')  -> '3a2bc'
    #        run_length_encode('annb$aa') -> 'a2nb$2a'
    #        run_length_encode('')        -> ''
    #        run_length_encode(BWT(dna_seq)) -> 151 chars become 143
```
9. Run-length decode
```
def run_length_decode(encoded):
    #1. Start with an empty output, and a place to hold digits not used yet
    decoded = ""
    digits = ""

    #2. Go through the encoded string one character at a time
    # a digit => part of a count, but more digits might follow, so hold onto it
    # a letter => the count in front of it is complete (none means 1): repeat the letter, clear the digits
    for each char in encoded:
        if char is a digit:
            digits = digits + char
        else:
            if digits is empty:
                repeat = 1
            else:
                repeat = int(digits)
            decoded = decoded + char repeated repeat times
            digits = ""

    #3. Return the rebuilt string
    return decoded

    # Check: run_length_decode('3a2bc')   -> 'aaabbc'
    #        run_length_decode('a2nb$2a') -> 'annb$aa'
    #        run_length_decode('12a')     -> 12 a's
    #        run_length_decode(run_length_encode(BWT(dna_seq))) == BWT(dna_seq) -> True
```

## Putting it together
10. Main
```
def main(reference, queries):
    #1. Show the reference
    print the reference and its length

    #2. Pattern search
    # find_match builds the suffix array, BWT, count and occur inside
    for each q in queries:
        print q and find_match(q, reference)

    #3. Compress the BWT and check the round trip
    text = reference + '$'
    bwt = BWT_from_suffix_array(text, suffix_array(text))
    encoded = run_length_encode(bwt)
    print the length before and after, and the encoded string
    print whether run_length_decode(encoded) == bwt

    # Check: main(generate_dna(150, seed=42), [reference[10:20], 'GAC', 'CATT', 'TTTT'])
    #        -> [10], [0, 57, 105, 123, 147], [69, 125, 130, 135], []; 151 -> 143 chars; True
```
## Checking the inverse (no separate inverse function)
```
Reversibility follows from the LF mapping. Starting in the $ row and repeatedly applying LF(row) = count[c] + occur[c][row] rebuilds the text backward, one character at a time.
The notebook does not run this walk, but every piece it needs is already in the code:
count[c] (where c’s block starts) comes from cal_count(), function 4
occur[c][row] (which occurrence c is) comes from cal_occur(), function 5
the same LF formula is used in update_range(), function 6, steps #3 and #4 (“LF of the top” and “LF of the bottom”), where it moves a whole block of rows instead of one row.
```
# Successes
Our biggest success was making sure the program is correct, not just that it runs. Instead of relying only on the starter examples, we checked our work in ways that would expose mistakes: two completely different methods of building the BWT produce the same result, the search handles every case the assignment asks for (including overlapping matches and patterns that are absent), and compressed data decodes back exactly. Using a reproducible DNA sequence means anyone can rerun the notebook and get the same results we report.\
We also came away understanding why the algorithm works. The last-to-first mapping turned out to be the idea that ties everything together, and once we understood it, both searching and reversing the transform made sense. Along the way, we caught problems the starter materials didn’t warn about, like a boundary case that breaks silently in Python, and we made the code reject bad input instead of producing wrong answers.

# Struggles

**Note on overall limitations of our program:**
- **Case-sensitive:** `'ANA'` does not find `'ana'`. Uppercase both first if case should not matter. Real genomes may be soft-masked (repeats in lowercase).
- **`cal_count()` does not skip `$`**, so it must be called on the reference without `$`.
- **`find_match()` adds its own `$`:** the reference must not contain `$`. `find_match('ana', 'banana$')` raises an error (the reference would end in two `$`).
- **Empty query:** `find_match('', reference)` returns every position, including the `$` position.
- **Index is rebuilt per query:** `find_match()` builds the suffix array, BWT, `count` and `occur` on every call.
- **Memory:** the naive BWT stores n rotations of length n; our simple `suffix_array()` also creates every suffix while sorting, so it is faster than the naive BWT but not truly O(n log n). Fast suffix sorting needs special algorithms. `occur` stores a count for every row.
- **Run-length encoding:** a text containing digits would break decoding; and on short inputs, or text without long runs, the encoding barely helps. On the random 150-bp DNA, the BWT went from 151 to 143 characters, because random DNA has little repeated structure. The paper finds the BWT needs large blocks to compress well - "a few kilobytes" (Burrows & Wheeler, 1994).

**Things we still don't understand:**

# Personal Reflections

## Group Leader

## Other members
Mara: The walkthrough provided by Marcus was extremely helpful in understanding the basis of the BWT. The principles behind the suffix array, doing the transform, and encoding/decoding were fairly straightforward once I sat down with a pencil and paper and did them by hand. The pattern-matching was the hardest part for me to understand conceptually. My group members were a fantastic team, and Trang was an excellent group leader. 

Trivedi: I started every function by working it out by hand on banana$ before writing any code, and that habit helped a lot. The BWT, suffix array, and count/occur tables made sense pretty quickly that way, but backward search took me the longest. I kept getting lost once the full sorted matrix was taken away, and it only clicked when I traced the LF mapping step by step and saw how count and occur jump from one range to the next. The part that surprised me most was how the same Python feature behaves differently across functions: text[-1] wrapping around was exactly what BWT_from_suffix_array needed, but occur[a][-1] silently broke update_range. Testing taught me a lot too. 'ana' is a palindrome, so it couldn't catch a search reading in the wrong direction, and building the occur table by hand is how I found the typo in the starter docstring. Comparing every result against a brute-force check on a random 150 bp sequence gave me real confidence the code was correct and not just running.

# Generative AI Appendix

**Generative AI:**\

Claude (Anthropic, Claude Opus 5.5) was used to debug run_length_decode and BWT_from_suffix_array, draft the test cases (including the random 150 bp DNA sequence), develop the main() function that brings everything together

**Citations**\
