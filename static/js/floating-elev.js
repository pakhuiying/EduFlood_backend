const formInputs = document.querySelectorAll(
  ".floating-contact-form .form-container .form-input .floating-elev .elev-container"
);

const contactIcon = document.querySelector(
  ".floating-contact-form .contact-icon .floating-elev .elev-icon"
);

const formContainer = document.querySelector(
  ".floating-contact-form .form-container .floating-elev .elev-container"
);

contactIcon.addEventListener("click", () => {
  formContainer.classList.toggle("active");
});

formInputs.forEach((i) => {
  i.addEventListener("focus", () => {
    i.previousElementSibling.classList.add("active");
  });
});

formInputs.forEach((i) => {
  i.addEventListener("blur", () => {
    if (i.value === "") {
      i.previousElementSibling.classList.remove("active");
    }
  });
});