
from ._dnasequencebank import DNASequenceBank

from ._models import (
        MaskedKmerModel,
        InterNucleotideDistanceModel,
        CodingDensityModel,
        Contig2VecModel,
        SequenceCoverageModel
    )

from ._sentence2vec.sentence2vec import Word, Sentence, get_word_frequency, sentence_to_vec

__version__ = 'dev.0.1'
__all__ = [
        'DNASequenceBank',
        'MaskedKmerModel',
        'InterNucleotideDistanceModel',
        'CodingDensityModel',
        'Contig2VecModel',
        'SequenceCoverageModel'
    ]
