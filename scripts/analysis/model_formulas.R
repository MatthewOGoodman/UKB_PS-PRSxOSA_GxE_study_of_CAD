formula_null_string_long =  paste(
  "~", 
  paste(
    model_vars_long,
    collapse = " + ") )

main_SA = paste(
  Exposure_Var
)


formula_null_string_spline =  paste(
  "~", 
  paste(
    model_vars_spline,
    collapse = " + ") )

formula_null_string = formula_null_string_spline

main_SA = paste(
  Exposure_Var
)

main_SA_fct = paste(
  Exposure_Var, "fct", sep = "_"
)


#Note: main_SA_spline is non-identifiable when SA is binary
# main_SA_spline = paste0( 
#   's(', main_SA,', df = 5)' 
# )

main_PGRS = paste(
  pathway
)

main_PGRS_spline = paste0( 
  's(', pathway,', df = 5)' 
)

interaction_SA_PGRS = paste(
  
  Exposure_Var,
  pathway,
  # paste0('I(', main_SA, ':', main_PGRS, ')'),
  paste0( main_SA, ':', main_PGRS ),
  
  sep = " + "
)

interaction_SA_PGRS_spline = paste(
  
  main_SA,
  # main_SA_spline,
  
  # main_PGRS_spline,
  
  # paste0('I(', main_SA_spline, ':', main_PGRS_spline, ')'),
  # paste0( main_SA_spline, ':', main_PGRS_spline ),
  # paste0( 's(', main_PGRS, ', by = ', main_SA, ', df = 5)' ),
  # paste0( 's( I(', main_SA, ' : ', main_PGRS, '), df = 5)' ),
  
  paste0( 's(', main_PGRS, ', by = ', main_SA_fct, ', k = 5)' ),
  # paste0( 's( I(', main_SA, ' : ', main_PGRS, '), k = 5)' ),
  
  sep = " + "
)


formula_SA_vs_PGRS = as.formula( paste( Exposure_Var, 
                                        
                                        paste(
                                          formula_null_string_long,
                                          main_PGRS,
                                          
                                          sep = " + "
                                        )
))

formula_SA_vs_PGRS_spline = as.formula( paste( Exposure_Var, 
                                               
                                               paste(
                                                 formula_null_string_spline,
                                                 main_PGRS_spline,
                                                 
                                                 sep = " + "
                                               )
))

formula_CAD_vs_PGRS = as.formula( paste( Outcome_Var, 
                                         
                                         paste(
                                           formula_null_string_long,
                                           main_PGRS,
                                           
                                           sep = " + "
                                         )
))

formula_CAD_vs_PGRS_spline = as.formula( paste( Outcome_Var, 
                                                
                                                paste(
                                                  formula_null_string_spline,
                                                  main_PGRS_spline,
                                                  
                                                  sep = " + "
                                                )
))


formula_CAD_vs_SA.PGRS = as.formula( paste( Outcome_Var, 
                                            
                                            paste(
                                              formula_null_string_long,
                                              interaction_SA_PGRS,
                                              
                                              sep = " + "
                                            )
))

formula_CAD_vs_SA.PGRS_spline = as.formula( paste( Outcome_Var, 
                                                   
                                                   paste(
                                                     formula_null_string_spline,
                                                     # interaction_SA_PGRS_spline,
                                                     interaction_SA_PGRS,
                                                     
                                                     sep = " + "
                                                   )
))

formula_CAD_vs_main.SA_main.PGRS = as.formula( paste( Outcome_Var, 
                                                             
                                                             paste(
                                                               formula_null_string_long,
                                                               main_SA,
                                                               main_PGRS,
                                                               
                                                               sep = " + "
                                                             )
))

formula_CAD_vs_main.SA_main.PGRS_spline = as.formula( paste( Outcome_Var, 
                                                             
                                                             paste(
                                                               formula_null_string_spline,
                                                               main_SA,
                                                               main_PGRS,
                                                               
                                                               sep = " + "
                                                             )
))

formula_CAD_vs_main.SA = as.formula( paste( Outcome_Var, 
                                                   
                                                   paste(
                                                     formula_null_string_long,
                                                     main_SA,

                                                     sep = " + "
                                                   )
))

formula_CAD_vs_main.SA_spline = as.formula( paste( Outcome_Var, 
                                                   
                                                   paste(
                                                     formula_null_string_spline,
                                                     main_SA,

                                                     sep = " + "
                                                   )
))
