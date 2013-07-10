# coding: utf-8
lib = File.expand_path('../lib', __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require 'FastFourier/fast-fourier-0.0.1.rb'

Gem::Specification.new do |spec|
  spec.name          = "FastFourier"
  spec.version       = FastFourier::VERSION
  spec.authors       = ["Justin W Smith"]
  spec.email         = ["justin.w.smith@gmail.com"]
  spec.description   = %q{An implementation of the Cooley-Tukey algorithm for the Discrete Fourier Transform}
  spec.summary       = %q{TODO: Write a gem summary}
  spec.homepage      = ""
  spec.license       = "MIT"

  spec.files         = `git ls-files`.split($/)
  spec.executables   = spec.files.grep(%r{^bin/}) { |f| File.basename(f) }
  spec.test_files    = spec.files.grep(%r{^(test|spec|features)/})
  spec.require_paths = ["lib"]

  spec.add_development_dependency "bundler", "~> 1.3"
  spec.add_development_dependency "rake"
end
