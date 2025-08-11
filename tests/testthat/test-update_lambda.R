# test_that("function works", {
#
#   G = 5
#   M = 3
#   n_basis = 10
#   n_id = 30
#   include_id_effect = FALSE
#   order = 1
#
#   D = diff(diag(n_basis), differences = order)
#   S = crossprod(D)
#
#   model_data = list(
#     G = G,
#     M = M,
#     n_basis = n_basis,
#     notpen_index = gen_notpen_index(n_basis = n_basis, M = M, order = order),
#     order = order,
#     n_id = n_id,
#     include_id_effect = FALSE,
#     D = D,
#     newD = D %*% t(D),
#     S = S,
#     basis_type = "bs"
#   )
#
#   alpha = matrix(rnorm(n_basis * M * (G - 1), mean = 0, sd = 1),
#                  nrow = M * n_basis)
#   alpha = cbind(alpha, 0)
#
#   new_lambda = update_lambda(
#     alpha = alpha,
#     a = 1,
#     b = 1,
#     model_data = model_data
#   ) %>%
#     expect_no_error() %>%
#     expect_no_message() %>%
#     expect_no_warning()
#
#   expect_equal(dim(new_lambda), c(M, G - 1))
#
#   new_lambda = update_lambda(
#     alpha = alpha,
#     a = 1,
#     b = 1,
#     model_data = model_data,
#     grouped_lambda = TRUE
#   ) %>%
#     expect_no_error() %>%
#     expect_no_message() %>%
#     expect_no_warning()
#
#   expect_equal(new_lambda[, 1], new_lambda[, 2])
#   expect_equal(new_lambda[, 1], new_lambda[, 3])
#
#   #### testing ps ####
#   model_data$basis_type = "ps"
#   new_lambda = update_lambda(
#     alpha = alpha,
#     a = 1,
#     b = 1,
#     model_data = model_data
#   ) %>%
#     expect_no_error() %>%
#     expect_no_message() %>%
#     expect_no_warning()
#
#
#
#
#   G = 5
#   M = 3
#   n_basis = 10
#   n_id = 30
#
#   model_data = list(
#     G = G,
#     M = M,
#     n_basis = n_basis,
#     notpen_index = gen_notpen_index(n_basis = n_basis, M = M),
#     n_id = n_id,
#     basis_type = "bs"
#   )
#
#   alpha_id = matrix(rnorm(n_id * (G - 1), mean = 0, sd = 1),
#                     nrow = n_id)
#   alpha_id = cbind(alpha_id, 0)
#
#   new_lambda = update_lambda_id(
#     alpha = alpha,
#     a = 1,
#     b = 1,
#     model_data = model_data
#   ) %>%
#     expect_no_error() %>%
#     expect_no_message() %>%
#     expect_no_warning()
#
#
#
# })
