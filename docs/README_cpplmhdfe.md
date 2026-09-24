# cpplmhdfe

## Validation and failure behavior

Fixed-effect and IRLS iteration limits and tolerances must be positive. A nonconverged projection or exhausted IRLS fit returns error 430. Joint Wald statistics use the retained regressor indices, so moving omitted columns does not change the reported test.
