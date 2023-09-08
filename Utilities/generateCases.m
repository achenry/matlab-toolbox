function [case_list, case_name_list, n_cases] = generateCases(case_basis, namebase, fast_fmt)
    
    % if fast_fmt is true, generate list of structs with fields-to-edit InflowWind,
    % Fst etc, each of which is a struct with fields-to-edit correspoding
    % to FAST parameters

    input_fields = fields(case_basis);

    if fast_fmt
        param_fields = struct;
        for f = 1:length(input_fields)
            sub_fields = fields(case_basis.(input_fields{f}));
            for ff = 1:length(sub_fields)
                if isstruct(case_basis.(input_fields{f}).(sub_fields{ff}))
                    subsub_fields = fields(case_basis.(input_fields{f}).(sub_fields{ff}));
                    param_fields.(input_fields{f}).(sub_fields{ff}) = subsub_fields;
                else 
                    param_fields.(input_fields{f}){ff} = sub_fields{ff};
                end
            end
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
            sub_fields = fields(case_basis.(input_fields{f}));

%             for p = 1:n_param_fields
            for ff = 1:length(sub_fields)
                if isstruct(case_basis.(input_fields{f}).(sub_fields{ff}))
                    param_basis = case_basis.(input_fields{f}).(sub_fields{ff});
                    n_param_fields = length(param_fields.(input_fields{f}).(sub_fields{ff}));
                    n_total_fields = n_total_fields + n_param_fields;

                    n = length(param_basis.(param_fields.(input_fields{f}).(sub_fields{ff})));
                    n_cases = n_cases * n;
                    n_vals = n_vals + n;
                    n_field_vals(f, ff) = n;
                    
                    numBits = n;
                    powersOf2 = 2.^(1-numBits:0);
                    bits = [];
                    for bit = 2.^(0:n-1)
                        bits = [bits; rem(floor(bit * powersOf2), 2)];
                    end
                    all_combinations = [repmat(all_combinations, length(bits), 1), ...
                        repmat(bits, max(size(all_combinations, 1), 1), 1)];
                else
                    param_basis = case_basis.(input_fields{f}).(sub_fields{ff});
                    n_param_fields = 1;%length(param_fields.(input_fields{f}));
                    n_total_fields = n_total_fields + n_param_fields;

%                     n = length(param_basis.(param_fields.(input_fields{f}){ff}));
                    n = length(param_basis);
                    n_cases = n_cases * n;
                    n_vals = n_vals + n;
                    n_field_vals(f, ff) = n;
                    
                    numBits = n;
                    powersOf2 = 2.^(1-numBits:0);
                    bits = [];
                    % new_combinations = zeros(max(size(all_combinations, 1), 1) * numBits, numBits);
                    clear new_combinations;
                    % new_combinations = [];
                    for bit = 2.^(0:n-1)
                        new_bits = rem(floor(bit * powersOf2), 2);
                        bits = [bits; new_bits];
                        bit_idx = size(bits, 1);
                        new_combination_start_idx = max(size(all_combinations, 1), 1) * (bit_idx - 1) + 1; % size(new_combinations, 1)
                        new_combination_end_idx = max(size(all_combinations, 1), 1) * bit_idx;
                        % repeat this new bit for each row in all_combinations
                        new_combinations(new_combination_start_idx:new_combination_end_idx, 1:size(new_bits, 2))...
                            = repmat(new_bits, max(size(all_combinations, 1), 1), 1);
                    end
                    all_combinations = [repmat(all_combinations, size(bits, 1), 1), ...
                                new_combinations];
                end
            end
        else
            n = length(case_basis.(input_fields{f}));
            n_cases = n_cases * n;
            n_vals = n_vals + n;
            n_field_vals(f) = n;
            
            numBits = n;
            powersOf2 = 2.^(1-numBits:0);
            bits = [];
            % new_combinations = zeros(max(size(all_combinations, 1), 1) * numBits, numBits);
            clear new_combinations;
            % new_combinations = [];
            for bit = 2.^(0:n-1)
                new_bits = rem(floor(bit * powersOf2), 2);
                bits = [bits; new_bits];
                bit_idx = size(bits, 1);
                new_combination_start_idx = max(size(all_combinations, 1), 1) * (bit_idx - 1) + 1; % size(new_combinations, 1)
                new_combination_end_idx = max(size(all_combinations, 1), 1) * bit_idx;
                % repeat this new bit for each row in all_combinations
                new_combinations(new_combination_start_idx:new_combination_end_idx, 1:size(new_bits, 2))...
                    = repmat(new_bits, max(size(all_combinations, 1), 1), 1);
            end
            all_combinations = [repmat(all_combinations, size(bits, 1), 1), ...
                        new_combinations];
                    % all_combinations = [repmat(all_combinations, length(bits), 1), ...
                    %     repmat(bits, max(size(all_combinations, 1), 1), 1)];
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