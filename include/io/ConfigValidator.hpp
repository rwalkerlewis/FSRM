#ifndef FSRM_CONFIG_VALIDATOR_HPP
#define FSRM_CONFIG_VALIDATOR_HPP

#include <string>
#include <vector>

namespace FSRM {

class ConfigReader;

// Strict configuration validator. Walks an already-parsed
// ConfigReader and rejects unknown sections, unknown keys within
// known sections, and a curated set of deprecated key / section
// forms left over from pass-12 housekeeping. Strict mode is the
// default; opt out by setting [META] strict_validation = false in
// the config (a stderr warning is still emitted on the opt-out).
//
// The validator is consulted from ConfigReader::loadFile after the
// raw key/value parse completes. Failures are reported via the
// returned ValidationResult; ConfigReader::loadFile then signals
// the failure to the caller (returns false in strict mode).
class ConfigValidator {
public:
    struct Result {
        bool valid = true;
        std::vector<std::string> errors;
        std::vector<std::string> warnings;
    };

    // Runs the strict validation pass. source_filename is woven
    // into error messages as the offending file path. Returns a
    // Result whose `valid` field is true when no errors were
    // accumulated. Warnings do not flip `valid`.
    static Result validate(const ConfigReader& reader,
                           const std::string& source_filename = "");

    // Convenience entry point that prints every error and warning
    // to stderr and returns Result::valid.
    static bool validateAndReport(const ConfigReader& reader,
                                  const std::string& source_filename = "");
};

} // namespace FSRM

#endif // FSRM_CONFIG_VALIDATOR_HPP
