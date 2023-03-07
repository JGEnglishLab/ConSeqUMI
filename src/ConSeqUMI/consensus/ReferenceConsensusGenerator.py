import re
from collections import Counter


class ReferenceConsensusGenerator:
    def __init__(self, *args, **kwargs):
        self.bufferLength = kwargs.get("bufferLength", 20)
        self.sequenceWindowLength = kwargs.get("sequenceWindowLength", 10)

    def justify_left_all_string_lengths_with_buffer(self, readSequences):
        maxLength = len(max(readSequences, key=len))
        justifiedSequences = [
            readSequence.ljust(maxLength + self.bufferLength)
            for readSequence in readSequences
        ]
        return justifiedSequences

    def initialize_consensus_sequence_front(self, readSequences):
        readSequenceFronts = [
            readSequence[: self.sequenceWindowLength] for readSequence in readSequences
        ]
        consensusFront = max(set(readSequenceFronts), key=readSequenceFronts.count)
        return consensusFront

    def find_next_character_in_sequence(self, readSequences, precursorSequence):
        precursorSequencePattern = re.compile(precursorSequence + ".")
        allCharacters = []
        for readSequence in readSequences:
            precursorSearchResult = precursorSequencePattern.search(readSequence)
            if not precursorSearchResult:
                continue
            readSequenceNextCharacter = precursorSearchResult.group(0)[-1]
            allCharacters.append(readSequenceNextCharacter)
        c = Counter(allCharacters)
        nextCharacter = max(set(allCharacters), key=allCharacters.count)
        return nextCharacter

    def generate_consensus_sequence(self, readSequences):
        justifiedReadSequences = self.justify_left_all_string_lengths_with_buffer(
            readSequences
        )
        consensusSequence = self.initialize_consensus_sequence_front(
            justifiedReadSequences
        )
        readSequenceLength = len(justifiedReadSequences[0])
        for i in range(readSequenceLength - self.bufferLength):
            readSequenceBufferWindows = [
                justifiedReadSequence[i : i + self.bufferLength]
                for justifiedReadSequence in justifiedReadSequences
            ]
            consensusSequenceWindow = consensusSequence[-self.sequenceWindowLength :]
            nextCharacter = self.find_next_character_in_sequence(
                readSequenceBufferWindows, consensusSequenceWindow
            )
            consensusSequence += nextCharacter
        return consensusSequence.strip(" ")
