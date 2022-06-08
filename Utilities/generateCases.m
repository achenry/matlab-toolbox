function [case_list, case_name_list, n_cases] = generateCases(case_basis, namebase, fast_fmt)
    
    % if fast_fmt is true, generate list of structs with fields-to-edit InflowWind,
    % Fst etc, each of which is a struct with fields-to-edit correspoding
    % to FAST parameters

    input_fields = fields(case_basis);

    if fast_fmt
        param_fields = struct;
        for f = 1:length(input_fields)
            param_fields.(input_fields{f}) = fields(case_basis.(input_fields{f}));
        end
    end

    n_cases = 1;
    n_vals = 0;
    n_total_fields = 0;
    n_field_vals = [];
    all_combinations = [];

    for f = 1:length(input_fields)
        n_field_vals = [n_field_vals; []];
        if fast_fmt
            param_basis = case_basis.(input_fields{f});
            n_param_fields = length(param_fields.(input_fields{f}));
            n_total_fields = n_total_fields + n_param_fields;
            for p = 1:n_param_fields
                n = length(param_basis.(param_fields.(input_fields{f}){p}));
                n_cases = n_cases * n;
                n_vals = n_vals + n;
                n_field_vals(f, p) = n;
                
                numBits = n;
                powersOf2 = 2.^(1-numBits:0);
                bits = [];
                for bit = 2.^(0:n-1)
                    bits = [bits; rem(floor(bit * powersOf2), 2)];
                end
                all_combinations = [repmat(all_combinations, length(bits), 1), ...
                    repmat(bits, max(size(all_combinations, 1), 1), 1)];
            end
        else
            n = length(case_basis.(input_fields{f}));
            n_cases = n_cases * n;
            n_vals = n_vals + n;
            n_field_vals(f) = n;
            
            numBits = n;
            powersOf2 = 2.^(1-numBits:0);
            bits = [];
            for bit = 2.^(0:n-1)
                bits = [bits; rem(floor(bit * powersOf2), 2)];
            end
            all_combinations = [repmat(all_combinations, length(bits), 1), ...
                    repmat(bits, max(size(all_combinations, 1), 1), 1)];
        end
    end
    
    % flatten cases into single array of cases
    case_list = [];
    case_name_list = {};
     
    for c = 1:n_cases
        new_struct = struct;

        if fast_fmt
            for f = 1:length(input_fields)
                for p = 1:length(param_fields.(input_fields{f}))
                    new_struct.(input_fields{f}).(param_fields.(input_fields{f}){p}) = 0;
                end
            end
        else
             for f = 1:length(input_fields)
                 new_struct.(input_fields{f}) = 0;

             end
        end
    
         case_list = [case_list, new_struct];
         case_name_list{c} = [namebase '_' num2str(c)];
    
    end
    

    for c = 1:n_cases
        bits = all_combinations(c, :);
        bin_i0 = 1;
        is_good_comb = true;
        for f = 1:length(input_fields)
            if fast_fmt
                for p = 1:length(param_fields.(input_fields{f}))
                    
                    bin_iend = bin_i0 + n_field_vals(f, p) - 1;
                    sub_bin_c = bits(bin_i0:bin_iend);
                    bin_i0 = bin_iend + 1;
                end
            else
                
                bin_iend = bin_i0 + n_field_vals(f) - 1;
                sub_bin_c = bits(bin_i0:bin_iend);
                bin_i0 = bin_iend + 1;
            end

            
        end
        
        bin_i0 = 1;
        for f = 1:length(input_fields)
            if fast_fmt
                for p = 1:length(param_fields.(input_fields{f}))
                    field_vals = case_basis.(input_fields{f}).(param_fields.(input_fields{f}){p});
                    bin_iend = bin_i0 + n_field_vals(f, p) - 1;
                    
                    sub_bin_c = bits(bin_i0:bin_iend);
                    if isnumeric(field_vals)
                        new_field_val = field_vals(sub_bin_c == 1);
                    elseif iscell(field_vals)
                        new_field_val = field_vals{sub_bin_c == 1};
                    end
                    
                    case_list(c).(input_fields{f}).(param_fields.(input_fields{f}){p}) = new_field_val;
                    bin_i0 = bin_iend + 1;
                end
            else
                field_vals = case_basis.(input_fields{f});
                bin_iend = bin_i0 + n_field_vals(f) - 1;
                sub_bin_c = bits(bin_i0:bin_iend);
                if isnumeric(field_vals)
                    new_field_val = field_vals(sub_bin_c == 1);
                elseif iscell(field_vals)
                    new_field_val = field_vals{sub_bin_c == 1};
                end
                case_list(c).(input_fields{f}) = new_field_val;
                bin_i0 = bin_iend + 1;
            end
        end
    end
end