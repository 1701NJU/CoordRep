# CoordRep Tools: legacy donor-field recovery diagnostics and syntax repair
#
# Torch-dependent modules are imported lazily so that pure-Python
# utilities (baselines, masking, ablation eval) work without GPU.

try:
    from .infer import mlm_topk, load_model_and_tokenizer
    from .tool_a_donor import DonorPredictor, evaluate_donor_prediction
    from .tool_b_repair import CoordRepRepairer, repair_coordrep
except ImportError:
    mlm_topk = None
    load_model_and_tokenizer = None
    DonorPredictor = None
    evaluate_donor_prediction = None
    CoordRepRepairer = None
    repair_coordrep = None

try:
    from .validate import is_valid_coordrep, validate_brackets, validate_structure
except ImportError:
    is_valid_coordrep = None
    validate_brackets = None
    validate_structure = None

from .baselines import RandomBaseline, FrequencyBaseline, ConditionalFrequencyBaseline

__all__ = [
    'mlm_topk',
    'load_model_and_tokenizer',
    'is_valid_coordrep',
    'validate_brackets',
    'validate_structure',
    'DonorPredictor',
    'evaluate_donor_prediction',
    'CoordRepRepairer',
    'repair_coordrep',
    'RandomBaseline',
    'FrequencyBaseline',
    'ConditionalFrequencyBaseline',
]
