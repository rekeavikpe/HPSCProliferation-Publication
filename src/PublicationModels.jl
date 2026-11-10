module PublicationModels
export publication_model, storage_model, model_label, publication_order
# Storage IDs are retained for compatibility with saved fits and checkpoints.
const PUBLICATION_BY_STORAGE = (0, 1, 2, 3, 4, 6, 5)
function publication_model(id::Integer)
    0 <= id <= 6 || throw(ArgumentError("Unknown storage model ID: $id"))
    PUBLICATION_BY_STORAGE[id + 1]
end
storage_model(id::Integer) = publication_model(id) # The permutation is its own inverse.
model_label(id::Integer) = "Model $(publication_model(id))"
publication_order() = sort(collect(0:6); by=publication_model)
end
