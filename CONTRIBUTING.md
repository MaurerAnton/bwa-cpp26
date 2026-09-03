# Contributing to BWA-CPP26

## Code Style

- C++26 standard
- No manual memory management (use Arena/Vector/PmrString)
- No exceptions in hot paths (use std::expected)
- const-correctness
- RAII for all resources

## Testing

```bash
# Run self-tests
./bwa-cpp26 test

# Run integration test
bash tests/integration/test_pipeline.sh
```

## Adding New Features

1. Write tests first (in main.cpp or new test file)
2. Implement feature
3. Run all tests
4. Update documentation

## Commit Messages

- `feat: description` - New feature
- `fix: description` - Bug fix
- `perf: description` - Performance improvement
- `refactor: description` - Code restructuring
- `docs: description` - Documentation
- `test: description` - Tests

## License

All contributions must be AGPL-3.0-only compatible.
