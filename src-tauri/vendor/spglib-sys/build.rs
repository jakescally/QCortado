use bindgen;
use cmake;

use std::env;
use std::path::PathBuf;

fn main() {
    println!("cargo:rerun-if-changed=wrapper.h");
    let dst = cmake::Config::new("spglib").build();
    for path in [
        dst.join("build"),
        dst.join("build").join("Release"),
        dst.join("build").join("RelWithDebInfo"),
        dst.join("build").join("MinSizeRel"),
        dst.join("lib"),
        dst.join("lib64"),
    ] {
        if path.exists() {
            println!("cargo:rustc-link-search=native={}", path.display());
        }
    }
    println!("cargo:rustc-link-lib=static=symspg");
    let bindings = bindgen::Builder::default()
        .header("wrapper.h")
        .parse_callbacks(Box::new(bindgen::CargoCallbacks))
        .generate()
        .expect("Unable to generate bindings");
    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Unable to write bindings");
}
