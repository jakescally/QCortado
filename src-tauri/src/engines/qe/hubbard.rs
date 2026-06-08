use super::types::{HubbardLrtResult, HubbardLrtValue};

fn normalize_float(raw: &str) -> Option<f64> {
    raw.replace('D', "E").replace('d', "e").parse::<f64>().ok()
}

fn is_element_symbol(raw: &str) -> bool {
    let mut chars = raw.chars();
    let Some(first) = chars.next() else {
        return false;
    };
    if !first.is_ascii_uppercase() {
        return false;
    }
    let rest: String = chars.collect();
    rest.len() <= 1 && rest.chars().all(|ch| ch.is_ascii_lowercase())
}

fn is_hubbard_manifold(raw: &str) -> bool {
    let mut chars = raw.chars().peekable();
    let mut saw_digit = false;
    while matches!(chars.peek(), Some(ch) if ch.is_ascii_digit()) {
        saw_digit = true;
        chars.next();
    }
    let Some(shell) = chars.next() else {
        return false;
    };
    saw_digit && matches!(shell, 's' | 'p' | 'd' | 'f') && chars.next().is_none()
}

fn split_target(raw: &str) -> Option<(String, String, String)> {
    let trimmed = raw
        .trim()
        .trim_matches(':')
        .trim_matches(',')
        .trim_start_matches("U(")
        .trim_end_matches(')');
    let mut pieces = trimmed.split('-');
    let element = pieces.next()?.trim();
    let manifold = pieces.next()?.trim();
    if pieces.next().is_some() || !is_element_symbol(element) || !is_hubbard_manifold(manifold) {
        return None;
    }
    Some((
        element.to_string(),
        manifold.to_string(),
        format!("{}-{}", element, manifold),
    ))
}

fn push_unique(values: &mut Vec<HubbardLrtValue>, value: HubbardLrtValue) {
    if let Some(existing) = values.iter_mut().find(|entry| entry.target == value.target) {
        *existing = value;
    } else {
        values.push(value);
    }
}

pub fn parse_hubbard_lrt_values(parameters_output: &str, hp_output: &str) -> Vec<HubbardLrtValue> {
    let mut values = Vec::new();
    parse_hubbard_parameter_table(parameters_output, &mut values);

    let combined = if values.is_empty() {
        format!("{}\n{}", parameters_output, hp_output)
    } else {
        parameters_output.to_string()
    };

    for line in combined.lines() {
        let normalized = line.trim();
        if normalized.is_empty() {
            continue;
        }
        let lower = normalized.to_ascii_lowercase();
        if !lower.contains('u') && !lower.contains("hubbard") {
            continue;
        }

        let tokens: Vec<&str> = normalized
            .split(|ch: char| ch.is_whitespace() || ch == '=')
            .filter(|token| !token.trim().is_empty())
            .collect();

        for (idx, token) in tokens.iter().enumerate() {
            let Some((element, manifold, target)) = split_target(token) else {
                continue;
            };
            let mut value_ev = None;
            for candidate in tokens.iter().skip(idx + 1) {
                let stripped = candidate.trim_matches(|ch: char| ch == ',' || ch == ';');
                if let Some(value) = normalize_float(stripped) {
                    value_ev = Some(value);
                    break;
                }
            }
            if let Some(value_ev) = value_ev {
                push_unique(
                    &mut values,
                    HubbardLrtValue {
                        element,
                        manifold,
                        target,
                        value_ev,
                    },
                );
            }
        }
    }

    values
}

fn parse_hubbard_parameter_table(text: &str, values: &mut Vec<HubbardLrtValue>) {
    let mut in_table = false;
    for line in text.lines() {
        let normalized = line.trim();
        if normalized.is_empty() {
            continue;
        }
        let lower = normalized.to_ascii_lowercase();
        if lower.contains("manifold") && lower.contains("hubbard u") {
            in_table = true;
            continue;
        }
        if !in_table || normalized.starts_with('=') || normalized.starts_with('-') {
            continue;
        }

        let tokens: Vec<&str> = normalized.split_whitespace().collect();
        if tokens.len() < 8 {
            continue;
        }
        let label = tokens[2];
        let manifold = tokens[6];
        let Some(value_ev) = tokens.last().and_then(|value| normalize_float(value)) else {
            continue;
        };
        if !is_element_symbol(label) || !is_hubbard_manifold(manifold) {
            continue;
        }
        push_unique(
            values,
            HubbardLrtValue {
                element: label.to_string(),
                manifold: manifold.to_string(),
                target: format!("{}-{}", label, manifold),
                value_ev,
            },
        );
    }
}

pub fn build_hubbard_lrt_result(
    hp_output: String,
    parameters_output: Option<String>,
    q_mesh: [u32; 3],
    artifacts: Vec<String>,
) -> HubbardLrtResult {
    let parameters_text = parameters_output.as_deref().unwrap_or("");
    let u_values = parse_hubbard_lrt_values(parameters_text, &hp_output);
    let converged = hp_output.contains("JOB DONE") || !u_values.is_empty();
    HubbardLrtResult {
        converged,
        q_mesh,
        u_values,
        raw_output: hp_output,
        parameters_output,
        artifacts,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_hubbard_targets_from_parameters_output() {
        let values = parse_hubbard_lrt_values(
            "
            Hubbard U parameters:
            U Fe-3d 5.4321
            U Ni-3d = 6.2500
            ",
            "",
        );
        assert_eq!(values.len(), 2);
        assert_eq!(values[0].target, "Fe-3d");
        assert!((values[0].value_ev - 5.4321).abs() < 1e-10);
        assert_eq!(values[1].target, "Ni-3d");
    }

    #[test]
    fn falls_back_to_hp_output() {
        let values = parse_hubbard_lrt_values("", "Hubbard U for Co-3d = 4.75 eV");
        assert_eq!(values[0].element, "Co");
        assert_eq!(values[0].manifold, "3d");
        assert_eq!(values[0].value_ev, 4.75);
    }

    #[test]
    fn parses_qe_hubbard_parameters_table() {
        let values = parse_hubbard_lrt_values(
            "
                                 Hubbard U parameters:

       site n.  type  label  spin  new_type  new_label  manifold  Hubbard U (eV)
         1        1   Ho       1      1         Ho         4f       6.1970
            ",
            "",
        );
        assert_eq!(values.len(), 1);
        assert_eq!(values[0].target, "Ho-4f");
        assert!((values[0].value_ev - 6.1970).abs() < 1e-10);
    }

    #[test]
    fn ignores_non_hubbard_hyphenated_words() {
        let values = parse_hubbard_lrt_values(
            "",
            "
            unit-cell volume = 300.5850
            kinetic-energy cutoff = 40
            U(Ho-4f) = 7.0000
            ",
        );
        assert_eq!(values.len(), 1);
        assert_eq!(values[0].target, "Ho-4f");
    }
}
