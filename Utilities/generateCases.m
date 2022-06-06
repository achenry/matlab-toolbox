function [case_suites_flat, n_cases] = generateCases(case_basis)

    case_fields = fields(case_basis);

    n_cases = 1;
    n_vals = 0;
    for f = 1:length(case_fields)
        n_cases = n_cases * length(getfield(case_basis, case_fields{f}));
        n_vals = n_vals + length(getfield(case_basis, case_fields{f}));
    end
    
    % flatten cases into single array of cases
    case_suites_flat = [];
    for c = 1:n_cases
        new_struct = struct;
    
         for f = 1:length(case_fields)
             new_struct.(case_fields{f}) = 0;%setfield(case_suites_flat(c), case_fields{f}, []);
         end
    
         case_suites_flat = [case_suites_flat, new_struct];
    
    end
    
    min_bin_digits = length(dec2bin(n_vals));
    combinations = zeros(n_cases, n_vals);
    c = 1;
    for v = 0:2^n_vals-1
        bin_char = dec2bin(v, n_vals);
        bin_arr = zeros(1, n_vals);
    
        for dig_idx = 1:n_vals
            bin_arr(dig_idx) = str2num(bin_char(dig_idx));
        end
    
        bin_arr = flip(bin_arr);
    
        bin_i0 = 1;
        is_good_comb = true;
        for f = 1:length(case_fields)
            n_field_vals = length(case_basis.(case_fields{f}));
            bin_iend = bin_i0 + n_field_vals - 1;
            sub_bin_c = bin_arr(bin_i0:bin_iend);
            bin_i0 = bin_iend + 1;
    
            if sum(sub_bin_c) ~= 1
                is_good_comb = false;
                break;
            end
        end
    
        if is_good_comb
            combinations(c, :) = bin_arr;
            bin_i0 = 1;
            for f = 1:length(case_fields)
                n_field_vals = length(case_basis.(case_fields{f}));
                bin_iend = bin_i0 + n_field_vals - 1;
                sub_bin_c = bin_arr(bin_i0:bin_iend);
    
                field_vals = case_basis.(case_fields{f});
                new_field_val = field_vals(find(sub_bin_c));
    
                % current_field_vals = getfield(case_suites_flat(c), case_fields{f});
    
                case_suites_flat(c).(case_fields{f}) = new_field_val;
    
                bin_i0 = bin_iend + 1;
            end
    
            c = c + 1;
        end
    end

end