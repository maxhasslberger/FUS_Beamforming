function [elementAll_pos, el2mask_ids, mask2el_ids] = el2mask_indexing(elementAll_pos_orig, n_arr_elements)

elementAll_pos = flip(elementAll_pos_orig, 1);

[elementAll_pos, el2mask_ids] = sortrows(elementAll_pos'); % refer element indices to mask -> right order in getDistributedSourceSignal

elementAll_pos = elementAll_pos';
elementAll_pos = flip(elementAll_pos, 1);

[~, mask2el_ids] = sort(el2mask_ids); % Refer back to original order
mask2el_ids = reshape(mask2el_ids, n_arr_elements, []); % Original order per transducer assuming each transducer has the same amount of elements

end
